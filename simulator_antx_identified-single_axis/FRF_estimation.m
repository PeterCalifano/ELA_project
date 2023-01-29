Nsamples = length(time_grid);

%% Compute and remove mean
q_mean = mean(q); % [rad/s]
ax_mean = mean(ax); % [m/s^2]
delta_mean = mean(Excit_signal);

% q_zm = q - q_mean;
% ax_zm = ax - ax_mean;
% delta_zm = Excit_signal - delta_mean;

%% Correlation function estimation

[R_qq, q_zm] = AutoCorrEst(q, Nsamples);
[R_aa, ax_zm] = AutoCorrEst(ax, Nsamples);
[R_dd, delta_zm] = AutoCorrEst(Excit_signal, Nsamples);


% R_xx = zeros(Nsamples, 1);
% 
% for nt = 1:Nsamples-1
%     for n = 1:Nsamples - nt
%         R_xx(nt) = R_xx(nt) + (x(n)*x(n + nt))/(Nsamples-nt);
%     end
% end

