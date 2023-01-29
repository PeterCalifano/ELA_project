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


