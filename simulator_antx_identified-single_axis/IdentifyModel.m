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
        Window = hanning(floor(N/K));

        [cohr1, f1] = mscohere(delta_zm, q_zm, Window, [], [], 1/sample_time); % f1 is in Hz
        [cohr2, ~] = mscohere(delta_zm, ax_zm, Window, [], [], 1/sample_time);

        figure
        subplot(2, 1, 1)
        semilogx(f1(cohr1 > 0.6), cohr1(cohr1 > 0.6), '-', 'LineWidth', 1.05);
        xlabel('Frequency [Hz]')
        grid minor
        title('Coherence $\delta$ --> $q$')

        subplot(2, 1, 2)
        semilogx(f1(cohr2 > 0.6), cohr2(cohr2 > 0.6), '-', 'LineWidth', 1.05);
        xlabel('Frequency [Hz]')
        grid minor
        title('Coherence $\delta$ --> $a_x$')

        f_window = f1(cohr1 > 0.6 & cohr2 > 0.6);
        % f_window = f_window(f_window <= 11);

        % Define data structure for greyest
        data_time = iddata([q_zm, ax_zm], delta_zm, sample_time);
        % Transform to frequency domain
        data = fft(data_time); % frequency in rad/s
        
        q_f = data.OutputData(:, 1);
        ax_f = data.OutputData(:, 2);
        delta_f = data.InputData;
        faxis = data.Frequency ./ (2*pi);

        figure;
        semilogy(faxis, abs(delta_f), 'k-', 'LineWidth', 1.01)
        xlabel('Frequency [Hz]')
        ylabel('$|M(f)|$ [dB]')
        grid minor
        axis auto
        title('Input signal - Frequency spectrum')

        f_lb = f_window(1);
        f_ub = f_window(end);

        id_f = faxis >= f_lb & faxis <= f_ub;
        q_data = q_f(id_f);
        ax_data = ax_f(id_f);
        delta_data = delta_f(id_f);

        data_to_fit = iddata([q_data, ax_data], delta_data, sample_time,...
            'Frequency', faxis(id_f), 'FrequencyUnit', 'Hz');

        % Guess parameters vector
        theta0 = th_true.*[1.2 1.1 0.8 1.1 0.7 1.1];

        % Generate initial guess for greyest function
        model_fun = 'LongDyn_ODE';
        [fitmodel, theta, est_unc] = greyest_wrapper(data_to_fit, model_fun, theta0);

        
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
                [fitmodel, theta, est_unc] = greyest_wrapper(data_to_fit, model_fun, theta0);

            case 1

                % Guess Model Transfer function
                TF_new = minreal(Hmodelstruct(theta0));

                % Determine frequency responses of H over grid f
                yH_sim = evalFreqR(TF_new, faxis, 'Hz');

                Nf = length(faxis);
                Nfcn = size(TF_new, 1)*size(TF_new, 2);

                % R = identity by default
                [J, eH] = J_LS(yH, yH_sim);
                eH = reshape(eH, Nf, 2*Nfcn);

                % Optimization procedure to get optimal theta parameters - Newton-Raphson

                CONVERGED = 0;
                Nparams = length(theta0);

                % Iterative search cycle
                c = 0; % Iteration counter
                Nmax = 15; % Max number of iterations
                DTHETA_THR = 1e-3; % Step tolerance
                DJ_THR = 1e-3; % Function Tolerance

                % Initialize guess vector
                theta = theta0;

                while ~CONVERGED

                    % Initialization of Gradient vector and Hessian Matrix
                    GJ = zeros(Nparams, 1);
                    HJ = zeros(Nparams, Nparams);

                    % Rinv must have size equal to 2*Number of entries in the transfer
                    % function matrix
                    Rinv = eye(2*Nfcn);

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
                    [J_new, e_new] = J_LS(yH, yH_sim_new);
                    % Reshape deviations array: Columns: Re1 Im1 Re2 Im2
                    e_new = reshape(e_new, Nf, 2*Nfcn);

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
                    fprintf('Iteration #%d\t cost function value %f ', c, J);
                    for ii = 1:Nparams
                        fprintf('\ttheta%d %f ', ii, theta(ii));
                    end
                    fprintf('\n');
                    c = c + 1;
                end


        end


end

% Error with respect to true theta
est_err = th_true' - theta;
fprintf("Initial Error: %4.4f\n", th_true - theta0);
fprintf("Error: %4.4f\n", est_err);


% %% CrossCorrelation fcn estimation
% tic
% [R_qd, R_dq] = CrossCorrEst(q_zm, delta_zm, Nsamples);
% [R_ad, R_da] = CrossCorrEst(ax_zm, delta_zm, Nsamples);
% [R_aq, R_qa] = CrossCorrEst(ax_zm, q_zm, Nsamples);
% toc

%% Signal split and windowing (overlapped)
% Reminder: trade-off to do with K number of intervals and M number of
% samples per interval --> bias vs variance of the estimate depending on
% the frequency bands of interest (control/determination?)

% Divide signal into K parts of length M from zero mean signals
% Windows Overlap coefficient
% x_frac = 0.1;
% % Number of Windows
% K = 10;
% 
% % Window length
% T_win = time_grid(end)/((K-1)*(1-x_frac)+1);
% 
% % Number of samples in each window
% N_win = round(Nsamples/((K-1)*(1-x_frac)+1));
% 
% % Subdivision od data into K records of individual lenght T_win
% 
% yax_int = cell(1,K);
% yq_int  = cell(1,K);
% xd_int  = cell(1,K);
% t_int   = cell(1,K);
% 
% yax_int{1} = ax_zm(1:N_win);
% yax_int{K} = ax_zm(end-N_win:end);
% 
% yq_int{1} = q_zm(1:N_win);
% yq_int{K} = q_zm(end-N_win:end);
% 
% xd_int{1} = delta_zm(1:N_win);
% xd_int{K} = delta_zm(end-N_win:end);
% t_int{1} = time_grid(1:N_win);
% t_int{K} = time_grid(end-N_win:end);
% 
% for k=2:K-1
% 
%     indexes = ceil((1-x_frac)*(k-1)*N_win:(1-x_frac)*(k-1)*N_win+N_win);
% 
%     while indexes(end) > length(time_grid)
%         indexes(end) = [];
%     end
% 
%     yax_int{k} = ax_zm(indexes);
%     yq_int{k} = q_zm(indexes);
%     xd_int{k} = delta_zm(indexes);
%     t_int{k} = time_grid(indexes);
% 
% end
% 
% % Create window functions and applied to each kth part
% 
% % Windowing
% ax_window = cell(1,K);
% q_window = cell(1,K);
% d_window = cell(1,K);
% 
% for k = 1:K
%     ax_window{k} = yax_int{k}.*bartlett(length(yax_int{k}));
%     q_window{k}  = yq_int{k}.*bartlett(length(yq_int{k}));
%     d_window{k}  = xd_int{k}.*bartlett(length(xd_int{k}));
% 
% end
% 
% %% DFT of signals
% % Apply DFT to each kth part
% 
% Y_ax = cell(1,K);
% Y_q  = cell(1,K);
% X_d  = cell(1,K);
% 
% for k = 1:K
%     X1 = fft(d_window{k}, Nsamples);
%     upper_index = floor((Nsamples+1)/2);
%     X_d{k} = X1(1:upper_index);
% 
%     Y1_ax = fft(ax_window{k},Nsamples);
%     Y_ax{k} = Y1_ax(1:upper_index);
% 
%     Y1_q = fft(q_window{k},Nsamples);
%     Y_q{k} = Y1_q(1:upper_index);
% 
% end
% 
% % Frequency
% f_axis = ((0:(Nsamples-1)/2)*1/sample_time/Nsamples)';
% 
% %% PSD rough estimate
% % Apply rough estimator to each kth segment --> Power * 2/T;
% 
% G_aa_rough = cell(1,K);
% G_qq_rough = cell(1,K);
% G_dd_rough = cell(1,K);
% G_da_rough = cell(1,K);
% G_dq_rough = cell(1,K);
% 
% for k = 1:K
%     G_aa_rough{k} = 2/T_win*abs(Y_ax{k}).^2 ;
%     G_qq_rough{k} = 2/T_win*abs(Y_q{k}).^2;
%     G_dd_rough{k} = 2/T_win*abs(X_d{k}).^2;
%     G_da_rough{k} = 2/T_win*conj(X_d{k}).*Y_ax{k};
%     G_dq_rough{k} = 2/T_win*conj(X_d{k}).*Y_q{k};
% 
% end
% 
% %% PSD smooth estimate
% % Apply either smooth estimate (mean of each rough estimate) or smooth-iterative
% % procedure (see slide 27 of PracticeClass7)
% 
% G_dd_mat = cell2mat(G_dd_rough);
% G_dd = mean(G_dd_mat,2);
% 
% G_aa_mat = cell2mat(G_aa_rough);
% G_aa = mean(G_aa_mat,2);
% 
% G_qq_mat = cell2mat(G_qq_rough);
% G_qq = mean(G_qq_mat,2);
% 
% G_da_mat = cell2mat(G_da_rough);
% G_da = mean(G_da_mat,2);
% 
% G_dq_mat = cell2mat(G_dq_rough);
% G_dq = mean(G_dq_mat, 2);
% 
% 
% % close all
% %
% % figure
% % loglog(f_axis, abs(G_dd))
% % grid minor
% %
% % figure
% % loglog(f_axis, abs(G_aa))
% % grid minor
% %
% % figure
% % loglog(f_axis, abs(G_qq))
% % grid minor
% %
% % figure
% % loglog(f_axis, abs(G_da))
% % grid minor
% %
% % figure
% % loglog(f_axis, abs(G_dq))
% % grid minor
% %
% % close all
% 
% 
% %% FRF estimation
% H1 = G_dq./G_dd;
% H2 = G_da./G_dd;
% 
% % Test Coherence evaluation function
% gamma2_dq = EstimateCoherence(G_dq, G_dd, G_qq, 2*pi*f_axis);
% title('Coherence of $H_1$')
% gamma2_da = EstimateCoherence(G_da, G_dd, G_aa, 2*pi*f_axis);
% title('Coherence of $H_2$')
% 
% % Select frequency ranges to use for Identification model
% gamma2_thr = 0.6;
% % Create Bool Mask
% MaskFreq = gamma2_dq >= gamma2_thr & gamma2_da >= gamma2_thr;
% % Extrac useful frequency points
% faxis_masked = f_axis(MaskFreq);
% % Estimated FRF of the system
% H_est = [H1(MaskFreq), H2(MaskFreq)];


%% Uncertainty assessment
% if method == 1
%     th_nom = [fitmodel.A(1, 1), fitmodel.A(2, 1), fitmodel.A(1, 2),...
%         fitmodel.A(2, 2), fitmodel.B(1), fitmodel.B(2)];
%     Hnom = Hmodelstruct(th_nom);
% 
%     y_nom = evalFreqR(Hnom(1), faxis_masked, 'Hz');
%     M = EstimateFisher(y_nom, th_nom, 1, faxis_masked);
%     Cth = M^-1;
%     [eigvec, eigval] = eig(Cth);
% end

% Alternatives for uncertainty assessment
% 1) From greyest --> uncertainty in terms of relative std. dev. of the
% parameters

% Print uncertainties
cell_param = {'Xu', 'Mu', 'Xq', 'Mq', 'Xd', 'Md'};
fprintf("\n\n");
for idp = 1:length(est_unc)
    fprintf("%s: std. dev.: %4.3g; relative dev.: %4.3g \n",...
        cell_param{idp}, est_unc(idp), est_unc(idp)./theta(idp));
end

% 2) Bodeplots of TFs from identified model, data through spafdr and
% transfer function estimated "manually" via spectral analysis

% Evaluate TF from "analytical" model (a)
TF_model = Hmodelstruct(theta);

% TF obtained with MATLAB function with identified parameters (b)
[Ai, Bi, Ci, Di] = LongDyn_ODE(theta(1), theta(2), theta(3), theta(4),theta(5), theta(6));
SS_model = ss(Ai, Bi, Ci, Di);
TF_from_ss = tf(SS_model);
% TF estimated with Spectral Analysis (K-W Theorem) (c)
EstimateFRF;
whos H_est faxis_masked

% NOT SURE WHY THE FREQUENCY AXIS IS SO SMALL IN BAND

w_axis = 2*pi*faxis_masked; % Frequency axis in rad/s
% Transfer function FRF of the identified model
[mag_a, phase_a, w_out] = bode(TF_model);
% Evaluate FRF starting from measurement data with MATLAB function
TF_spafdr = spafdr(iddata([q_data, ax_data], delta_data, sample_time, 'Frequency', faxis(id_f), 'FrequencyUnit', 'Hz'));
[mag_b, phase_b] = bode(TF_spafdr, w_out);

% NOTE: 20 or 10 to convert H in log scale --> check what bode() uses
legend_cell = {'greyest model', 'spafdr', 'KW estimator', '', ''};
% TF: delta --> q
id = 1;
figure;
% Magnitude plots
subplot(2, 1, 1);
hold on;
semilogx(w_out, squeeze(mag_a(1, 1, :)), '-', 'LineWidth', 1.05, 'DisplayName', 'greyest model');
semilogx(w_out, squeeze(mag_b(1, 1, :)), '-', 'LineWidth', 1.05, 'DisplayName', 'spafdr');
semilogx(w_axis, 20*log10(abs(H_est(:, 1))), '-', 'LineWidth', 1.05, 'DisplayName', 'KW estimator');
grid minor
axis auto;

ax = gca;
ax.XAxisLocation = 'bottom'; 
ax.YAxisLocation = 'left';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.LineWidth = 1.04;
xlabel('Frequency [rad/s]')
ylabel('Gain [dB]')
title('Bode magnitude of $H_{\delta q}$')
legend()

hold off;

% Phase
subplot(2, 1, 2);
hold on;
semilogx(w_out, (squeeze(phase_a(id, 1, :))), '-', 'LineWidth', 1.05, 'DisplayName', 'greyest model');
semilogx(w_out, (squeeze(phase_b(id, 1, :))), '-', 'LineWidth', 1.05, 'DisplayName', 'spafdr');
semilogx(w_axis, rad2deg(angle(H_est(:, id))), '-', 'LineWidth', 1.05, 'DisplayName', 'KW estimator');
yline(-180, 'k--', 'DisplayName', '');
yline(+180, 'k--', 'DisplayName', '');
grid minor
axis auto;

ax = gca;
ax.XAxisLocation = 'bottom'; 
ax.YAxisLocation = 'left';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.LineWidth = 1.04;
xlabel('Frequency [rad/s]')
ylabel('Phase [deg]')
title('Bode phase of $H_{\delta q}$')
legend(legend_cell);


hold off;

% TF: delta --> ax
id = 2;
figure;
% Magnitude plots
subplot(2, 1, 1);
hold on;
semilogx(w_out, squeeze(mag_a(id, 1, :)), '-', 'LineWidth', 1.05, 'DisplayName', 'greyest model');
semilogx(w_out, squeeze(mag_b(id, 1, :)), '-', 'LineWidth', 1.05, 'DisplayName', 'spafdr');
semilogx(w_axis, 20*log10(abs(H_est(:, id))), '-', 'LineWidth', 1.05, 'DisplayName', 'KW estimator');
grid minor;
axis auto;

ax = gca;
ax.XAxisLocation = 'bottom'; 
ax.YAxisLocation = 'left';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.LineWidth = 1.04;
xlabel('Frequency [rad/s]');
ylabel('Gain [dB]');
title('Bode magnitude of $H_{\delta a_x}$')
legend();

hold off;

% Phase
subplot(2, 1, 2);
hold on;
semilogx(w_out, (squeeze(phase_a(id, 1, :))), '-', 'LineWidth', 1.05, 'DisplayName', 'greyest model');
semilogx(w_out, (squeeze(phase_b(id, 1, :))), '-', 'LineWidth', 1.05, 'DisplayName', 'spafdr');
semilogx(w_axis, rad2deg(angle(H_est(:, id))), '-', 'LineWidth', 1.05, 'DisplayName', 'KW estimator');
yline(-180, 'k--');
yline(+180, 'k--');
grid minor;
axis auto;

ax = gca;
ax.XAxisLocation = 'bottom'; 
ax.YAxisLocation = 'left';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.LineWidth = 1.04;
xlabel('Frequency [rad/s]');
ylabel('Phase [deg]');
title('Bode phase of $G_{\delta a_x}$');
legend(legend_cell);

hold off;


% 3) Bodeplots of identified model with confidence intervals from (1)
figure;
b_plot = bodeplot(fitmodel);
showConfidence(b_plot, 1);
grid minor
axis auto
ax = gca;
ax.XAxisLocation = 'bottom';
ax.YAxisLocation = 'left';
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.LineWidth = 1.04;
title('Bode diagrams $[G_{\delta q}; G_{\delta a_x}]$ with 1\sigma confidence interval', 'Interpreter', 'Latex')

% 3b) Pole-Zero with uncertainty interval
p = figure;
pz_plot = iopzplot(fitmodel);
showConfidence(pz_plot, 3)
hold on;

for i = 1:2
%     subplot(2, 1, i)
    hold on;
    grid minor
    axis auto
    ax = gca;
    ax.XAxisLocation = 'bottom';
    ax.YAxisLocation = 'left';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.LineWidth = 1.05;
end
title('Poles/Zeros of $H_{\delta a_x}$ with 3$\sigma$ confidence interval')

% 4) MCM to sample uncertainty space and get corresponding FRF
% Assume Gaussian distribution
N_samples = 50;
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

figure;
bodeplot(fitmodel, '*-r'); % , 'Linewidth', 1.04, 'Color', '#dd3322');
hold on;
bodeplot(unc_sys, '-k'); %, 'Linewidth', 1.02', 'Color', '#111111');
grid minor
axis auto
legend('Nominal identified model', 'Uncertain models')

figure;
iopzplot(fitmodel, '*-r'); % , 'Linewidth', 1.04, 'Color', '#dd3322');
hold on;
iopzplot(unc_sys, '-k'); %, 'Linewidth', 1.02', 'Color', '#111111');
grid minor
axis auto
legend('Nominal identified model', 'Uncertain models')

% For validation:
% 1) Identified model simulation with parameters and corresponding
% uncertainty --> get new simulated measurements --> re-identification
% (re-identified parameters must be similar for the process to be
% validated)
% 2) Simulated again with identified model and check consistency with
% reference responses (with the exception of the noise)

% Estimation of the variance q, ax from zero input portions of the output
% signals?







