%% Options

% 0: MATLAB greyest, 1: K-W for PSD --> optim_method
method = 1;

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
        
        K = 15;
        Window = hanning(floor(N/K));

        [cohr1, f1] = mscohere(delta_zm, q_zm, Window, [], [], 1/sample_time);
        [cohr2, ~] = mscohere(delta_zm, ax_zm, Window, [], [], 1/sample_time);

        figure
        subplot(2, 1, 1)
        semilogx(f1(cohr1 > 0.6), cohr1(cohr1 > 0.6), '-');
        xlabel('Frequency [Hz]')
        grid minor

        subplot(2, 1, 2)
        semilogx(f1(cohr2 > 0.6), cohr2(cohr2 > 0.6), '-');
        xlabel('Frequency [Hz]')
        grid minor

        f_window = f1(cohr1 > 0.6 & cohr2 > 0.6);
        % f_window = f_window(f_window <= 11);

        % Define data structure for greyest
        data_time = iddata([q_zm, ax_zm], delta_zm, sample_time);
        % Transform to frequency domain
        data = fft(data_time);

        q_f = data.OutputData(:, 1);
        ax_f = data.OutputData(:, 2);
        delta_f = data.InputData;
        faxis = data.Frequency ./ 2*pi;

        f_lb = f_window(1);
        f_ub = f_window(end);

        id_f = faxis >= f_lb & faxis <= f_ub;
        q_data = q_f(id_f);
        ax_data = ax_f(id_f);
        delta_data = delta_f(id_f);

        data_to_fit = iddata([q_data, ax_data], delta_data, sample_time, 'Frequency', faxis(id_f));

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

return

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
x_frac = 0.1;
% Number of Windows
K = 10;

% Window length
T_win = time_grid(end)/((K-1)*(1-x_frac)+1);

% Number of samples in each window
N_win = round(Nsamples/((K-1)*(1-x_frac)+1));

% Subdivision od data into K records of individual lenght T_win

yax_int = cell(1,K);
yq_int  = cell(1,K);
xd_int  = cell(1,K);
t_int   = cell(1,K);

yax_int{1} = ax_zm(1:N_win);
yax_int{K} = ax_zm(end-N_win:end);

yq_int{1} = q_zm(1:N_win);
yq_int{K} = q_zm(end-N_win:end);

xd_int{1} = delta_zm(1:N_win);
xd_int{K} = delta_zm(end-N_win:end);
t_int{1} = time_grid(1:N_win);
t_int{K} = time_grid(end-N_win:end);

for k=2:K-1

    indexes = ceil((1-x_frac)*(k-1)*N_win:(1-x_frac)*(k-1)*N_win+N_win);

    while indexes(end) > length(time_grid)
        indexes(end) = [];
    end

    yax_int{k} = ax_zm(indexes);
    yq_int{k} = q_zm(indexes);
    xd_int{k} = delta_zm(indexes);
    t_int{k} = time_grid(indexes);

end

% Create window functions and applied to each kth part

% Windowing
ax_window = cell(1,K);
q_window = cell(1,K);
d_window = cell(1,K);

for k = 1:K
    ax_window{k} = yax_int{k}.*bartlett(length(yax_int{k}));
    q_window{k}  = yq_int{k}.*bartlett(length(yq_int{k}));
    d_window{k}  = xd_int{k}.*bartlett(length(xd_int{k}));

end

%% DFT of signals
% Apply DFT to each kth part

Y_ax = cell(1,K);
Y_q  = cell(1,K);
X_d  = cell(1,K);

for k = 1:K
    X1 = fft(d_window{k}, Nsamples);
    X_d{k} = X1(1:(Nsamples+1)/2);

    Y1_ax = fft(ax_window{k},Nsamples);
    Y_ax{k} = Y1_ax(1:(Nsamples+1)/2);

    Y1_q = fft(q_window{k},Nsamples);
    Y_q{k} = Y1_q(1:(Nsamples+1)/2);

end

% Frequency
f_axis = ((0:(Nsamples-1)/2)*1/sample_time/Nsamples)';

%% PSD rough estimate
% Apply rough estimator to each kth segment --> Power * 2/T;

G_aa_rough = cell(1,K);
G_qq_rough = cell(1,K);
G_dd_rough = cell(1,K);
G_da_rough = cell(1,K);
G_dq_rough = cell(1,K);

for k = 1:K
    G_aa_rough{k} = 2/T_win*abs(Y_ax{k}).^2 ;
    G_qq_rough{k} = 2/T_win*abs(Y_q{k}).^2;
    G_dd_rough{k} = 2/T_win*abs(X_d{k}).^2;
    G_da_rough{k} = 2/T_win*conj(X_d{k}).*Y_ax{k};
    G_dq_rough{k} = 2/T_win*conj(X_d{k}).*Y_q{k};

end

%% PSD smooth estimate
% Apply either smooth estimate (mean of each rough estimate) or smooth-iterative
% procedure (see slide 27 of PracticeClass7)

G_dd_mat = cell2mat(G_dd_rough);
G_dd = mean(G_dd_mat,2);

G_aa_mat = cell2mat(G_aa_rough);
G_aa = mean(G_aa_mat,2);

G_qq_mat = cell2mat(G_qq_rough);
G_qq = mean(G_qq_mat,2);

G_da_mat = cell2mat(G_da_rough);
G_da = mean(G_da_mat,2);

G_dq_mat = cell2mat(G_dq_rough);
G_dq = mean(G_dq_mat, 2);


% close all
%
% figure
% loglog(f_axis, abs(G_dd))
% grid minor
%
% figure
% loglog(f_axis, abs(G_aa))
% grid minor
%
% figure
% loglog(f_axis, abs(G_qq))
% grid minor
%
% figure
% loglog(f_axis, abs(G_da))
% grid minor
%
% figure
% loglog(f_axis, abs(G_dq))
% grid minor
%
% close all


%% FRF estimation
H1 = G_dq./G_dd;
H2 = G_da./G_dd;

% Test Coherence evaluation function
% gamma2_dq = EstimateCoherence(G_dq, G_dd, G_qq, 2*pi*f_axis);
% title('Coherence of $H_1$')
% gamma2_da = EstimateCoherence(G_da, G_dd, G_aa, 2*pi*f_axis);
% title('Coherence of $H_2$')

% Select frequency ranges to use for Identification model
gamma2_thr = 0.6;
% Create Bool Mask
MaskFreq = gamma2_dq >= gamma2_thr & gamma2_da >= gamma2_thr;
% Extrac useful frequency points
faxis_masked = f_axis(MaskFreq);
% Estimated FRF of the system
H_est = [H1(MaskFreq), H2(MaskFreq)];


return

%% Uncertainty assessment
if method == 1
    th_nom = [fitmodel.A(1, 1), fitmodel.A(2, 1), fitmodel.A(1, 2),...
        fitmodel.A(2, 2), fitmodel.B(1), fitmodel.B(2)];
    Hnom = Hmodelstruct(th_nom);

    y_nom = evalFreqR(Hnom(1), faxis_masked, 'Hz');
    M = EstimateFisher(y_nom, th_nom, 1, faxis_masked);
    Cth = M^-1;
    [eigvec, eigval] = eig(Cth);
end



