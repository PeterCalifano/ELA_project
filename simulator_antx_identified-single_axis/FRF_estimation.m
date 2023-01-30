Nsamples = length(time_grid);

%% Compute and remove mean
% q_mean = mean(q); % [rad/s]
% ax_mean = mean(ax); % [m/s^2]
% delta_mean = mean(Excit_signal);

%% AutoCorrelation fcn estimation
[R_qq, q_zm, q_mean] = AutoCorrEst(q, Nsamples);
[R_aa, ax_zm, ax_mean] = AutoCorrEst(ax, Nsamples);
[R_dd, delta_zm, delta_mean] = AutoCorrEst(Excit_signal, Nsamples);

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
x_frac = 0.5;
K = 111;

% Window length
T_win = time_grid(end)/((K-1)*(1-x_frac)+1);

% Number of samples in each window
N_win = round(Nsamples/((K-1)*(1-x_frac)+1));

% Subdivision od data into K records of individual lenght T_win

yax_int = cell(1,K);
yq_int  = cell(1,K);
xd_int  = cell(1,K);
t_int   = cell(1,K);

for k=2:K-1
    yax_int{1} = ax_zm(1:N_win);
    yax_int{k} = ax_zm((1-x_frac)*(k-1)*N_win:(1-x_frac)*(k-1)*N_win+N_win);
    yax_int{K} = ax_zm(end-N_win:end);

    yq_int{1} = q_zm(1:N_win);
    yq_int{k} = q_zm((1-x_frac)*(k-1)*N_win:(1-x_frac)*(k-1)*N_win+N_win);
    yq_int{K} = q_zm(end-N_win:end);

    xd_int{1} = delta_zm(1:N_win);
    xd_int{k} = delta_zm((1-x_frac)*(k-1)*N_win:(1-x_frac)*(k-1)*N_win+N_win);
    xd_int{K} = delta_zm(end-N_win:end);

    t_int{1} = time_grid(1:N_win);
    t_int{k} = time_grid((1-x_frac)*(k-1)*N_win:(1-x_frac)*(k-1)*N_win+N_win);
    t_int{K} = time_grid(end-N_win:end);
    
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
    X1 = fft(d_window{k},Nsamples);
    X_d{k} = X1(1:(Nsamples+1)/2);

    Y1_ax = fft(ax_window{k},Nsamples);
    Y_ax{k} = Y1_ax(1:(Nsamples+1)/2);

    Y1_q = fft(q_window{k},Nsamples);
    Y_q{k} = Y1_q(1:(Nsamples+1)/2);

end

% frequency
f = ((0:(Nsamples-1)/2)*1/sample_time/Nsamples)';

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
G_dq = mean(G_dq_mat,2);
%%
close all
figure
loglog(f,abs(G_dd))
grid minor

figure
loglog(f,abs(G_aa))
grid minor

figure
loglog(f,abs(G_qq))
grid minor


figure
loglog(f,abs(G_da))
grid minor

figure
loglog(f,abs(G_dq))
grid minor




%% FRF estimation


