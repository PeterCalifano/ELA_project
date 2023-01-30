Nsamples = length(time_grid);

%% Compute and remove mean
% q_mean = mean(q); % [rad/s]
% ax_mean = mean(ax); % [m/s^2]
% delta_mean = mean(Excit_signal);

%% AutoCorrelation fcn estimation
[R_qq, q_zm, q_mean] = AutoCorrEst(q, Nsamples);
[R_aa, ax_zm, ax_mean] = AutoCorrEst(ax, Nsamples);
[R_dd, delta_zm, delta_mean] = AutoCorrEst(Excit_signal, Nsamples);

%% CrossCorrelation fcn estimation
tic
[R_qd, R_dq] = CrossCorrEst(q_zm, delta_zm, Nsamples);
[R_ad, R_da] = CrossCorrEst(ax_zm, delta_zm, Nsamples);
[R_aq, R_qa] = CrossCorrEst(ax_zm, q_zm, Nsamples);
toc

%% Signal split and windowing (overlapped)
% Reminder: trade-off to do with K number of intervals and M number of
% samples per interval --> bias vs variance of the estimate depending on
% the frequency bands of interest (control/determination?)

% Divide signal into K parts of length M from zero mean signals
x_frac = 0.5;
K = 111;

% Window length
T_win = time_grid/((K-1)*(1-x_frac)+1);

% Number of samples in each window
N_win = round(N/((K-1)*(1-x_frac)+1));

% Subdivision od data into K records of individual lenght T_win

x_int = cell(1,K);
return
y_int

% Create window functions and applied to each kth part


%% DFT of signals
% Apply DFT to each kth part


%% PSD rough estimate
% Apply rough estimator to each kth segment --> Power * 2/T;

%% PSD smooth estimate
% Apply either smooth estimate (mean of each rough estimate) or smooth-iterative
% procedure (see slide 27 of PracticeClass7)


%% FRF estimation


