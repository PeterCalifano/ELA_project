%% Options

% 0: MATLAB greyest, 1: K-W for PSD --> optim_method
method = 0;

% 0: MATLAB greyest or 1: Newton-Raphson
optim_method = 0;

Nsamples = length(time_grid);

%% AutoCorrelation fcn estimation
[R_qq, q_zm, q_mean] = AutoCorrEst(q, Nsamples);
[R_aa, ax_zm, ax_mean] = AutoCorrEst(ax, Nsamples);
[R_dd, delta_zm, delta_mean] = AutoCorrEst(Mtot, Nsamples);

switch method

    case 0 % MATLAB Built-in Greyest identification

        %% Compute and remove mean
        N = length(delta_zm);
        
        K = 8;
        % Create window for averaging
        Window = hanning(floor(N/K));

        % Cross-correlation of input and output signals in time-domain
        [cohr1, f1] = mscohere(delta_zm, q_zm, Window, [], [], 1/sample_time); % f1 is in Hz
        [cohr2, ~] = mscohere(delta_zm, ax_zm, Window, [], [], 1/sample_time);


        % Coherence plot (OK)
        time_coherence = figure;
        semilogx(f1, cohr1, '-', 'LineWidth', 1.03);
        hold on;
        semilogx(f1, cohr2, '-', 'LineWidth', 1.03);
        f_window = f1(cohr1 > 0.6 & cohr2 > 0.6);
        
        f_lb = f_window(1);
        f_window = f_window(f_window <= 15);
        f_ub = f_window(end);

        xline(f_lb, 'k--', num2str(f_lb) + " Hz", 'LineWidth', 1.04, 'FontSize', 12, ...
            'LabelOrientation', 'horizontal' , 'LabelVerticalAlignment','bottom');

        xline(f_ub, 'k--', num2str(f_ub) + " Hz", 'LineWidth', 1.04, 'FontSize', 12, ...
            'LabelOrientation', 'horizontal', 'LabelVerticalAlignment','bottom');

        title('Coherence - Output and Input time series')
        xlabel('Frequency [Hz]')
        ylabel('$|\gamma_{uy}|^2$')
        legend('$\delta$ with $q$', '$\delta$ with $a_x$', '', '', 'Location', 'bestoutside');

        % Default Options
        grid minor
        axis auto
        ylim([0 1.1])
        ax_gca = gca;
        ax_gca.XAxisLocation = 'bottom';
        ax_gca.YAxisLocation = 'left';
        ax_gca.XMinorTick = 'on';
        ax_gca.YMinorTick = 'on';
        ax_gca.LineWidth = 1.04;
        hold off;

        % Define data structure for greyest
        data_time = iddata([q_zm, ax_zm], delta_zm, sample_time, ...
            'InputName', {'M_{tot}'}, 'OutputName', {'q','ax'});
        % Transform to frequency domain
        data_freq = fft(data_time); % frequency in rad/s

        q_f = data_freq.OutputData(:, 1);
        ax_f = data_freq.OutputData(:, 2);
        delta_f = data_freq.InputData;
        % f_axis in [Hz] Validated by looking at Excitation signal 
        faxis = data_freq.Frequency./(2*pi); 

        % Input spectrum plot (OK)
        figure;
        loglog(faxis, abs(delta_f), 'k-', 'LineWidth', 1.02)
        xlabel('Frequency [Hz]')
        ylabel('$|M(\omega)|$ [dB]')
        title('Input signal $M_{tot}$ - Frequency spectrum');
        % Default Options
        grid minor
        axis auto
        ax_gca = gca;
        ax_gca.XAxisLocation = 'bottom';
        ax_gca.YAxisLocation = 'left';
        ax_gca.XMinorTick = 'on';
        ax_gca.YMinorTick = 'on';
        ax_gca.LineWidth = 1.04;

        % Extract samples at frequencies of interest
        id_f = faxis >= f_lb & faxis <= f_ub;

        q_data = q_f(id_f);
        ax_data = ax_f(id_f);
        delta_data = delta_f(id_f);

        % Encapsulate frequency domain data for fitting
        data_to_fit = iddata([q_data, ax_data], delta_data, sample_time,...
            'Frequency', faxis(id_f), 'FrequencyUnit', 'Hz',...
            'InputName', {'M_{tot}'}, 'OutputName', {'q','ax'}, ...
            'OutputUnit', {'rad/s', 'm/s2'});

        % Guess parameters vector
        theta0 = th_true.*[1.2 1.1 0.8 1.1 0.7 1.1];

        % Generate initial guess for greyest function
        model_fun = 'LongDyn_ODE';
        [model, theta, est_unc] = greyest_wrapper(data_to_fit, model_fun, theta0, 0);
        
    case 1

        % Call script to estimate FRF of the system by K-W Theory (PSD)
        EstimateFRF; % Output: H_est, faxis_masked
        whos H_est;

        % Determine data for identification with Newton-Raphson method
        faxis = faxis_masked; % [Hz]
        yH = formatFRF(H_est);

        % Guess parameters vector
        theta0 = th_true.*[1.2 1.1 0.8 1.1 0.7 1.1];

        switch optim_method
            case 0

                % Assembly Frequency Response Data from H
                freqdata = zeros(2, 1, length(faxis_masked));
                freqdata(1, 1, :) = H_est(:, 1);
                freqdata(2, 1, :) = H_est(:, 2);

                data_to_fit = frd(freqdata, faxis_masked, sample_time, 'FrequencyUnit', 'Hz');

                % Generate initial guess for greyest function
                model_fun = 'LongDyn_ODE';
                [model, theta, est_unc] = greyest_wrapper(data_to_fit, model_fun, theta0, 1);

            case 1

                % Guess Model Transfer function
                TF_new = minreal(Hmodelstruct(theta0));

                % Determine frequency responses of H over grid f
                yH_sim = evalFreqR(TF_new, faxis, 'Hz');

                Nf = length(faxis);
                Nfcn = size(TF_new, 1)*size(TF_new, 2);

                % R = identity by default
                [J, eH, R] = J_LS(yH, yH_sim);
                eH = reshape(eH, Nf, 2*Nfcn);

                % Optimization procedure to get optimal theta parameters - Newton-Raphson

                CONVERGED = 0;
                Nparams = length(theta0);

                % Iterative search cycle
                c = 0; % Iteration counter
                Nmax = 15; % Max number of iterations
                DTHETA_THR = 1e-9; % Step tolerance
                DJ_THR = 1e-12; % Function Tolerance

                % Initialize guess vector
                theta = theta0;
                % Rinv must have size equal to 2*Number of entries in the transfer
                % function matrix
                Rinv = R^(-1);

                while ~CONVERGED

                    % Initialization of Gradient vector and Hessian Matrix
                    GJ = zeros(Nparams, 1);
                    HJ = zeros(Nparams, Nparams);

                    % idp indexes the parameter in theta
                    dFRFdth = ComputeSensitivity(yH_sim, faxis, theta);

                    for ff = 1:Nf
                        dFRFdth_idp = zeros(4, Nparams);
                        for idp = 1:Nparams
                            dFRFdth_idp(:, idp) = dFRFdth{idp}(ff, :);
                        end

                        t = -eH(ff, :) * Rinv * dFRFdth_idp;
                        GJ = GJ + t';

                        HJ = HJ + dFRFdth_idp' * Rinv * dFRFdth_idp;
                    end

                    eigH = eig(HJ);

                    if min(eigH) <= 0
                        warning('Hessian not positive definite')
                    end

                    % Evaluate dtheta to find new guess vector
                    diff_theta = -HJ\GJ;

                    if ~iscolumn(theta)
                        theta = theta';
                    end

                    % Compute new guess vector
                    theta_new = theta + diff_theta;
                    % Evaluate new Transfer Function
                    TF_new = minreal(Hmodelstruct(theta_new));
                    % Determine frequency response
                    yH_sim_new = evalFreqR(TF_new, faxis, 'Hz');
                    % Evaluate new cost and deviations
                    [J_new, e_new, R] = J_LS(yH, yH_sim_new);
                    % Reshape deviations array: Columns: Re1 Im1 Re2 Im2
                    e_new = reshape(e_new, Nf, 2*Nfcn);

                    Rinv = R^(-1);

                    % Stop the iteration if the cost is increasing
                    cost_prev = J;
                    cost_now = J_new;
                    if cost_now - cost_prev > 1e3
                        error('Divergence occurred: Cost is increasing')
                    end

                    % check convergence
                    if (norm(diff_theta) < DTHETA_THR ...
                            || norm(J_new - J) < DJ_THR ...
                            || c > Nmax)
                        CONVERGED = true;
                    end

                    % update
                    theta = theta_new;
                

                    J = J_new;
                    eH = e_new;
                    yH_sim = yH_sim_new;

                    % DEBUG
                    fprintf('Iteration #%d cost function value %f ', c, J);
                    for ii = 1:Nparams
                        fprintf('\ttheta%d %f ', ii, theta(ii));
                    end
                    fprintf('\n');
                    c = c + 1;
                end
        end
end


%% Uncertainty assessment

% Error with respect to true theta
est_err = abs(th_true' - theta);
est_err0 = abs(th_true - theta0);

% Print uncertainties
fprintf("Initial guess error and fit model error:" + ...
    "\nXu: %3.5f; \t %3.5f" + ...
    "\nMu: %3.5f; \t %3.5f" + ...
    "\nXq: %3.5f; \t %3.5f" + ...
    "\nMq: %3.5f; \t %3.5f" + ...
    "\nXd: %3.5f; \t %3.5f" + ...
    "\nMd: %3.5f; \t %3.5f\n",...
    est_err0(1), est_err(1), ...
    est_err0(2), est_err(2), ...
    est_err0(3), est_err(3), ...
    est_err0(4), est_err(4), ...
    est_err0(5), est_err(5), ...
    est_err0(6), est_err(6));


fprintf("\nIdentified model parameters and relative 3-sigma uncertainty (Gaussian distr. assumption):" + ...
    "\nXu = %3.3g\t 3sig %%: %3.3g" + ...
    "\nMu = %3.3g\t 3sig %%: %3.3g" + ...
    "\nXq = %3.3g\t 3sig %%: %3.3g" + ...
    "\nMq = %3.3g\t 3sig %%: %3.3g" + ...
    "\nXd = %3.3g\t 3sig %%: %3.3g" + ...
    "\nMd = %3.3g\t 3sig %%: %3.3g\n\n", ...
    theta(1), 100*abs(3*est_unc(1)./theta(1)), ...
    theta(2), 100*abs(3*est_unc(2)./theta(2)), ...
    theta(3), 100*abs(3*est_unc(3)./theta(3)), ...
    theta(4), 100*abs(3*est_unc(4)./theta(4)), ...
    theta(5), 100*abs(3*est_unc(5)./theta(5)), ...
    theta(6), 100*abs(3*est_unc(6)./theta(6)));


%% Bodeplots of TFs (greyest, spafdr, estimator)
% TF obtained with MATLAB function with greyest identified parameters (a)
[Ai, Bi, Ci, Di] = LongDyn_ODE(theta(1), theta(2), theta(3), theta(4),theta(5), theta(6));
SS_model = ss(Ai, Bi, Ci, Di);
TF_model = tf(SS_model);
% TF estimated with spafdr from measurement data with MATLAB function
TF_spafdr = spafdr(iddata([q_data, ax_data], delta_data, sample_time, 'Frequency', 2*pi*faxis(id_f), 'FrequencyUnit', 'rad/s'));
% TF estimated with Spectral Analysis (K-W Theorem) (c)
run('EstimateFRF');
whos H_hat faxis_masked

% Frequency axis in rad/s from FRF estimators
w_axis = 2*pi*faxis_masked; 
% Transfer function FRF of the identified model
[mag_a, phase_a, w_out] = bode(TF_model, w_axis);
[mag_b, phase_b, w_out] = bode(TF_spafdr, w_out);

% Phase wrapping
phase_a = wrapTo180(phase_a);
phase_b = wrapTo180(phase_b);

legend_cell = {'greyest model', 'spafdr()', 'FRF estimator', '', ''};

% TF: delta --> q (1) and delta --> ax (2)
for id = 1:2
figure;
% Magnitude plots
subplot(2, 1, 1);
semilogx(w_out, squeeze(mag_a(id, 1, :)), '-', 'LineWidth', 1.05, 'DisplayName', 'greyest model');
hold on;
semilogx(w_out, squeeze(mag_b(id, 1, :)), '-', 'LineWidth', 1.05, 'DisplayName', 'spafdr()');
semilogx(w_axis, 20*log10(abs(H_hat(:, id))), '-', 'LineWidth', 1.05, 'DisplayName', 'FRF estimator');

xlabel('Frequency [rad/s]')
ylabel('Magnitude [dB]')

if id == 1
    title('Magnitude of $H_{\delta q}$')
elseif id == 2
    title('Magnitude of $H_{\delta a_x}$')
end

legend('Location', 'bestoutside')
% Default Options
grid minor
axis auto
ax_gca = gca;
ax_gca.XAxisLocation = 'bottom';
ax_gca.YAxisLocation = 'left';
ax_gca.XMinorTick = 'on';
ax_gca.YMinorTick = 'on';
ax_gca.LineWidth = 1.04;
hold off;

% Phase
subplot(2, 1, 2);
semilogx(w_out, (squeeze(phase_a(id, 1, :))), '-', 'LineWidth', 1.05, 'DisplayName', 'greyest model');
hold on;
semilogx(w_out, (squeeze(phase_b(id, 1, :))), '-', 'LineWidth', 1.05, 'DisplayName', 'spafdr');
semilogx(w_axis, rad2deg(angle(H_hat(:, id))), '-', 'LineWidth', 1.05, 'DisplayName', 'FRF estimator');
yline(-180, 'k--',  '-180°', 'FontSize', 12, 'LineWidth', 1.03, 'LabelHorizontalAlignment', 'Left');
yline(+180, 'k--',  '+180°', 'FontSize', 12, 'LineWidth', 1.03, 'LabelHorizontalAlignment', 'Left', 'LabelVerticalAlignment', 'bottom');

xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
if id == 1
    title('Phase of $H_{\delta q}$')
elseif id == 2
    title('Phase of $H_{\delta a_x}$');
end
legend(legend_cell, 'Location', 'bestoutside');
% Default Options
grid minor
axis auto
ax_gca = gca;
ax_gca.XAxisLocation = 'bottom';
ax_gca.YAxisLocation = 'left';
ax_gca.XMinorTick = 'on';
ax_gca.YMinorTick = 'on';
ax_gca.LineWidth = 1.04;
hold off;

end


%% Fitting comparison
% NOTE: compare uses FIT metric
compare_plot = figure;
Data = iddata([q_data, ax_data], delta_data, sample_time, 'Frequency', 2*pi*faxis(id_f),...
    'InputName', {'M_{tot}'}, 'OutputName', {'q','ax'});
compare(Data, model, 'r');
title('Evaluation of goodness of fit from greyest identified model', 'Interpreter', 'Latex');

% Default options applied to all axes handles
axes_handles = findall(compare_plot, 'type', 'axes');
set(axes_handles, 'YMinorTick', 'on', 'XMinorTick', 'on', 'LineWidth', 1.04, ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on');
hold off;
legend('Interpreter', 'Latex')

%% Bodeplots of identified model with confidence intervals from (1)
bode_fig = figure;
b_plot = bodeplot(model);
showConfidence(b_plot, 3);
grid minor
axis auto
title('Bode diagrams $[H_{\delta q}; H_{\delta a_x}]$ with 3$\sigma$ confidence interval',...
    'Interpreter', 'Latex')
% Default options applied to all axes handles
axes_handles = findall(bode_fig, 'type', 'axes');
set(axes_handles, 'YMinorTick', 'on', 'XMinorTick', 'on', 'LineWidth', 1.04, ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on');
hold off;

% 3b) Pole-Zero with uncertainty interval
pfig = figure;
pz_plot = iopzplot(model);
showConfidence(pz_plot, 3)
hold on;

title('Poles/Zeros of $[H_{\delta q}; H_{\delta a_x}]$ with 3$\sigma$ confidence interval',...
    'Interpreter', 'Latex');
% Default options applied to all axes handles
axes_handles = findall(pfig, 'type', 'axes');
set(axes_handles, 'YMinorTick', 'on', 'XMinorTick', 'on', 'LineWidth', 1.04, ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on');
hold off;


%% MCM to assess uncertainty space
% Assume Gaussian distribution
N_samples = 10;
C = diag(est_unc.^2);
% Sample uncertain parameter space
Theta_sampled = mvnrnd(theta, C, N_samples)';

% Build uncertain state space
A_unc = zeros(3, 3, N_samples);
B_unc = zeros(3, 1, N_samples);
C_unc = zeros(2, 3, N_samples);
D_unc  = zeros(2, 1, N_samples);

for id = 1:N_samples
    [A_unc(:, :, id), B_unc(:, :, id), C_unc(:, :, id), D_unc(:, :, id)] = LongDyn_ODE(Theta_sampled(1, id), Theta_sampled(2, id),Theta_sampled(3, id),Theta_sampled(4, id),...
        Theta_sampled(5, id), Theta_sampled(6, id));
end

unc_sys = ss(A_unc, B_unc, C_unc, D_unc);

% Define input/output names
unc_sys.u = 'M';
unc_sys.y = {'q','ax'};

waxis = 0:0.5:150;
bode_MCM = figure;
nominalbode = bodeplot(model, '.-b', waxis); 
showConfidence(nominalbode, 3);
hold on;
bodeplot(unc_sys, 'k-', waxis); 
title('MCM Bode diagrams of $[H_{\delta q}; H_{\delta a_x}]$',...
    'Interpreter', 'Latex');
% Default options applied to all axes handles
axes_handles = findall(bode_MCM, 'type', 'axes');
set(axes_handles, 'YMinorTick', 'on', 'XMinorTick', 'on', 'LineWidth', 1.04, ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on');
hold off;
legend('Nominal identified model', 'Uncertain models', 'Location', 'best')

pfig_MCM = figure;
nominalpz = iopzplot(model, 'r*'); 
% showConfidence(nominalpz, 3);
hold on;
iopzplot(unc_sys, 'k'); 
title('MCM Poles/Zeros of $[H_{\delta q}; H_{\delta a_x}]$',...
    'Interpreter', 'Latex');
% Default options applied to all axes handles
axes_handles = findall(pfig_MCM, 'type', 'axes');
set(axes_handles, 'YMinorTick', 'on', 'XMinorTick', 'on', 'LineWidth', 1.04, ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on');
hold off;
legend('Nominal identified model', 'Uncertain models', 'Location', 'best')






