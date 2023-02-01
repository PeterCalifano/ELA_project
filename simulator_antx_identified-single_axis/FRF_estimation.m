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


close all

figure
loglog(f_axis, abs(G_dd))
grid minor

figure
loglog(f_axis, abs(G_aa))
grid minor

figure
loglog(f_axis, abs(G_qq))
grid minor

figure
loglog(f_axis, abs(G_da))
grid minor

figure
loglog(f_axis, abs(G_dq))
grid minor

close all


%% FRF estimation
H1 = G_dq./G_dd;
H2 = G_da./G_dd;

% Test Coherence evaluation function
gamma2_dq = EstimateCoherence(G_dq, G_dd, G_qq, 2*pi*f_axis);
gamma2_da = EstimateCoherence(G_da, G_dd, G_aa, 2*pi*f_axis);

yH = formatFRF([H1, H2]);

% H1real = real(H1);
% H1imag = imag(H1);
% 
% H2real = real(H2);
% H2imag = imag(H2);

% "Measurements" vector
% yH = [[H1real; H1imag], [H2real; H2imag]];
 
% Determine TF model class
% First alternative: by looking at state space form of the system
% Parameters in theta: Xu Xq Mu Mq Xd Md 
% syms th [6, 1]
% syms s % TF variable
% 
% % A = [Xu Xq -9.81; Mu Mq 0; 0 1 0];
% % B = [Xd; Md; 0];
% % C = [0 1 0; Xu Xq 0];
% % D = [0;Xd];
% 
% A = [th(1) th(2) -9.81; th(3), th(4) 0; 0 1 0];
% B = [th(5); th(6); 0];
% C = [0, 1, 0; th(1), th(2), 0];
% D = [0; th(5)];
% 
% Phi = (s*eye(3) - A)^(-1);
% H = C*Phi*B+D;
% 
% % Display matrix
% pretty(H);

th_true = [-0.1068, 0.1192, -5.9755,  -2.6478, -10.1647, 450.71];

% Test function Hmodel
fcn_true = minreal(Hmodelstruct(th_true));

% %%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%
% figure;
% bode(G(1), f);
% hold on
% bode(fcn(1), f);
% grid minor
% legend;
% 
% figure;
% bode(G(2), f);
% hold on
% bode(fcn(2), f);
% grid minor
% legend;
%%%%%%%%%%%%%%%%%%%%%%%% DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%

TF = minreal(Hmodelstruct(th_true));
% Determine frequency responses of H over grid f
yH_sim = evalFreqR(TF, f_axis, 'Hz');

Nf = length(f_axis);

% Assume noise variance matrix (identity for now)
% 1st col: re(G1), 2nd col: im(G1), 3rd col: re(G2), 4th col: im(G2)
% eH = reshape(yH - yH_sim, length(f_axis), 4);

% J = 0;
% for id = 1:length(f_axis)
%     R = eye(4);
%     % Evaluate cost function
%     % J = 1/2 * sum(e'* (R^-1) * e);
%     J = J + 0.5 * eH(id, :) * R^-1 * eH(id, :)';
% end

[J, eH] = J_LS(yH, yH_sim);

dTFdth = ComputeSensitivity(yH_sim, f_axis, th_true);

% Optimization procedure to get optimal theta parameters
% Newton-Raphson 
% TO DO: 
% 2.1) Verify Sensitivity Function
% 2.2) Code function to execute Newton Raphson scheme
% 3) Test output

return

FLAG_CONVERGENCE = 0;

% iteration
while ~FLAG_CONVERGENCE

    % Compute gradient and hessian of J wrt theta
%     s = computeFRFsensitivity(omega_vec, G_th, theta);


    G = zeros(nth, 1);
    H = zeros(nth, nth);

    for kk=1:K
        dydtheta = zeros(2, nth);
        for ii = 1:nth
            dydtheta(:, ii) = s(ii).dy(:, kk);
        end
        t = -eH(:,kk)' * Rinv * dydtheta;
        G = G + t';

        H = H + dydtheta' * Rinv * dydtheta;
    end

    % compute step (Newton-Raphson)
    delta_theta = -H\G;

    theta_new = theta + delta_theta;

    [J_new, e_new] = costFunction(g_m, omega_vec, G_th, theta_new);

    % check convergence
    if (norm(delta_theta) < TRESHOLD_THETA ...
            || norm(J_new - J) < TRESHOLD_J ...
            || c > Nmax)
        FLAG_CONVERGENCE = true;
    end

    % update
    theta = theta_new;
    c = c+1;
    J = J_new;
    eH = e_new;

    % DEBUG
    fprintf('Iteration #%d\t cost function value %f ',c,J);
    for ii=1:nth
        fprintf('\ttheta%d %f ', ii, theta(ii));
    end
    fprintf('\n');

end


%% Function



