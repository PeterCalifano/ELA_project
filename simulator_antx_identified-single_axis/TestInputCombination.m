close all
clear
clc

%% Load workspace and parameters
addpath('common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');

DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize);

rng shuffle;

% NOTE: parameters required to the simulator and configuration are loaded
% and set inside the function

comb = 2;
switch comb
    case 1
        load('comb1_best.mat') % LSW + RBS
    case 2
        load('comb2_best.mat') % LSW + 3211
end

% Options
display_flag = 1;
metric_selector = 1; % sum(variances)

[J, fitmodel, estimates, sigma] = IdentificationExperiment(optimal_input, theta0, signal_type, metric_selector, display_flag);
FIT = fitmodel.Report.Fit.FitPercent;

%%
fprintf("\n-------- Identified model parameters and relative 3-sigma uncertainty (Gaussian distr. assumption) --------" + ...
    "\nXu = %3.3g\t 3sig%%: %3.3g" + ...
    "\nMu = %3.3g\t 3sig%%: %3.3g" + ...
    "\nXq = %3.3g\t 3sig%%: %3.3g" + ...
    "\nMq = %3.3g\t 3sig%%: %3.3g" + ...
    "\nXd = %3.3g\t 3sig%%: %3.3g" + ...
    "\nMd = %3.3g\t 3sig%%: %3.3g", ...
    estimates(1), 100*abs(3*sigma(1)./estimates(1)), ...
    estimates(2), 100*abs(3*sigma(2)./estimates(2)), ...
    estimates(3), 100*abs(3*sigma(3)./estimates(3)), ...
    estimates(4), 100*abs(3*sigma(4)./estimates(4)), ...
    estimates(5), 100*abs(3*sigma(5)./estimates(5)), ...
    estimates(6), 100*abs(3*sigma(6)./estimates(6)));

fprintf('\n\n');
fprintf('-------- Errors with respect to "truth" model --------');

bias_perc = 100 * abs((estimates' - th_true)./th_true);

cell_param = {'Xu', 'Mu', 'Xq', 'Mq', 'Xd', 'Md'};
fprintf("\n");
for idp = 1:length(sigma)
    fprintf("%s %% estimation error: %4.3g \n",...
        cell_param{idp}, bias_perc(idp));
end

% Generate Validation signal
rng(1);
params.dt = 0.001;
params.f0 = 0;
params.ff = 15;
params.t0 = 0;
params.tf = 150;

[signal, timevec] = GenerateInput(params, 1);
ExcitationM = [timevec, signal];
%
% Get reference output
sim_obj = SetModel(th_true, ExcitationM);
output = sim(sim_obj);

% Signals Pre-Processing
N_delay = 1;
[Mtot_ref, ax_ref, q_ref, time_grid_ref, pitch_angle_ref] = OutputPreProcess(output, N_delay);

clear output sim_obj

% Get identified model output
rng(1);
sim_obj = SetModel(estimates, ExcitationM);
output = sim(sim_obj);

% Signals Pre-Processing
N_delay = 1;
[Mtot, ax, q, time_grid, pitch_angle] = OutputPreProcess(output, N_delay);

figure;
plot(Mtot);
hold on;
plot(Mtot_ref);

figure;
plot(q);
hold on;
plot(q_ref);


% Time-domain signals errors
output_error(:, 1) = q_ref - q;
output_error(:, 2) = ax_ref - ax;
output_error(:, 3) = pitch_angle_ref - pitch_angle;

subplot(2, 1, 1)
plot(time_grid_ref, (output_error(:, 1)), '-', 'Linewidth', 1.02);
xlabel('Time [s]')
ylabel('Pitch rate [rad/s]')
title('Pitch rate errors - Time domain');
grid minor
axis auto;
legend();

subplot(2, 1, 2)
plot(time_grid_ref, (output_error(:, 2)), '-', 'Linewidth', 1.02);
title('Acceleration errors - Time domain');
grid minor
axis auto;
legend();


% Frequency Response with bode plots
figure;
w_axis = 0:0.05:(2*pi*100);
bodeplot(Hmodelstruct(th_true), 'r-',...
    Hmodelstruct(estimates), 'k-', w_axis);

legend('"Truth" model', 'Identified model');
grid minor;

% Pole-Zero plot
p = figure;
iopzplot(Hmodelstruct(th_true), 'r');
hold on;
pz_plot = iopzplot(fitmodel, 'k');
showConfidence(pz_plot, 3)
legend('"Truth" model', 'Identified model', '', '');
grid minor;

%% Computation of Validation Indexes
evalFIT = @(qmeas, q) 100 * max([0, 1 - (norm(qmeas - q)).^2 ./ norm(qmeas - mean(qmeas))]);
evalPEC = @(qmeas, q) 1/sqrt(length(qmeas)) .* (norm(qmeas - q)).^2;
evalVAF = @(qmeas, q) 100 * max([0, 1 - var(qmeas - q)./var(qmeas)]);

ModelFIT = [evalFIT(pitch_angle_ref, pitch_angle); evalFIT(q_ref, q); evalFIT(ax_ref, ax)];
ModelPEC = [evalPEC(pitch_angle_ref, pitch_angle); evalPEC(q_ref, q); evalPEC(ax_ref, ax)];
ModelVAF = [evalVAF(pitch_angle_ref, pitch_angle); evalVAF(q_ref, q); evalVAF(ax_ref, ax)];

% Display results:
fprintf("\n\nModelFIT: theta: %3.1g%%, q: %3.1g%%, ax: %3.1g%%\n  ", ModelFIT(1), ModelFIT(2), ModelFIT(3));
fprintf("ModelPEC: theta: %3.1g, q: %3.1g, ax: %3.1g\n  ", ModelPEC(1), ModelPEC(2), ModelPEC(3));
fprintf("ModelVAF: theta: %3.1g%%, q: %3.1g%%, ax: %3.1g%%\n  ", ModelVAF(1), ModelVAF(2), ModelVAF(3));

