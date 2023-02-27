%% Options
results_only = 1;
if ~exist("sample_time", "var")
    sample_time = 1/250; % [s]
end

%% Experiment optimization
theta0 = th_true.*[1.2 1.1 0.8 1.1 0.7 1.1];

comb = 2;
completed_flag = 0;

switch comb
    case 1 % Linear Sine Sweep + RBS
        % Pameters: 1) K, 2) tf, 3) f0, 4) ff, 5) T
        % Initial guess for optimization
        x0 = [8, 295, 0.05, 20, 1];

        LB = [5, 1*60, 0.001, 1, 0.05];
        UB = [40, 5*60, 0.2, 150, 10];

        signal_type = [1, 3];

        % Percentage of time for the first signal type
        params.tfrac = 0.5;

    case 2 % Linear Sine Sweep + 3211 (HALF)
        % Pameters: 1) K, 2) tf, 3) f0, 4) ff, 5) N
        % Initial guess for optimization
        x0 = [8, 200, 0.05, 10, 6];

        LB = [5, 1*60, 0, 0.8, 1];
        UB = [50, 5*60, 0.05, 50, 20];
        signal_type = [1, 4];

        % Percentage of time for the first signal type
        params.tfrac = 0.5;

    case 3 % Sine Sweep 
        % Pameters: 1) K, 2) tf, 3) f0, 4) ff
        % Initial guess for optimization
        x0 = [8, 295, 0.05, 20];

        LB = [5, 1*60, 0.001, 1];
        UB = [40, 5*60, 0.2, 150];

        signal_type = 1;

    case 4 % Linear Sine Sweep + 3211 (80/20)
        
        % Pameters: 1) K, 2) tf, 3) f0, 4) ff, 5) N
        % Initial guess for optimization
        x0 = [8, 295, 0.05, 20, 1];

        LB = [5, 1*60, 0.001, 1, 0.05];
        UB = [40, 5*60, 0.2, 150, 10];

        signal_type = [1, 3];
        % Percentage of time for the first signal type
        params.tfrac = 0.8; 

    case 5 % RBS
        % Pameters: 1) K, 2) tf, 3) T
        % Initial guess for optimization
        x0 = [8, 295, 1];

        LB = [5, 1*60, 0.005];
        UB = [40, 5*60, 10];

        signal_type = 3;

    case 6 % 3211
        % Pameters: 1) K, 2) tf, 3) N
        % Initial guess for optimization
        x0 = [8, 295, 5];

        LB = [5, 1*60, 2];
        UB = [40, 5*60, 30];

        signal_type = 4;

    case 7 % Doublet
        % Pameters: 1) K, 2) tf, 3) N
        % Initial guess for optimization
        x0 = [8, 295, 5];

        LB = [5, 1*60, 2];
        UB = [40, 5*60, 50];

        signal_type = 5;

end

if results_only == 0 
% x0 = [8, 295, 0.02, 50];
% x0 = [8,295,20];
metric_selector = 1;

J = IdentificationExperiment(x0, theta0, signal_type, metric_selector, 1);
disp("Initial J cost: " + num2str(J));


% Set-up options
fmincon_opts = optimoptions("fmincon", "Algorithm", 'interior-point',...
       'FunctionTolerance', 1e-5, 'MaxIterations', 15,...
       'UseParallel', false, 'Display', 'iter',...
       'OptimalityTolerance', 1e-10, 'ConstraintTolerance', 1e-10);

ga_opts = optimoptions("ga", "Display", 'iter', 'CrossoverFraction', 0.7, ...
    'FunctionTolerance', 1e-4, 'MaxTime', 6*3600, 'PopulationSize', 50,...
    'MaxStallGenerations', 15,'UseParallel', false);

% Run optimization
% [optimal_input, ] = fmincon(@(x) IdentificationExperiment(x, theta0,...
%     metric_selector, 0), x0, [], [], [], [], LB, UB, [], fmincon_opts);

[optimal_input] = ga(@(x) IdentificationExperiment(x, theta0, signal_type,...
    metric_selector, 0), length(LB), [], [], [], [], LB, UB, [], [], ga_opts);

[optimal_input2] = fmincon(@(x) IdentificationExperiment(x, theta0, signal_type,...
    metric_selector, 0), optimal_input, [], [], [], [], LB, UB, [], fmincon_opts);

save('comb2modified_best.mat');
completed_flag = 1;

end

%% Input test and validation
if completed_flag == 1

    % Load and Assign optimal input combinations
    [optimal_input, signal_type, params] = LoadOptimalInput(comb, sample_time, optimal_input);

elseif results_only == 1 && completed_flag == 0

    % Test all the optimal input combinations if only results are desired
    J_cell = nan(7, 1);
    fitmodel_cell = cell(7, 1);
    est_params = cell(7, 1);
    est_unc = cell(7, 1);
    dataset = cell(7, 1);
    optimal_input = cell(7, 1);
    signal_type = cell(7, 1);
    amplitudes_cell = cell(7, 1);
    excitation_cell = cell(7, 1);
    FIT_perc = cell(7, 1);

    est_err = nan(6, 7);

    name_cell = ["LSW + RBS", "LSW + 3211 (50-50)", "LSW", "LSW + 3211 (80-20)", "RBS", "3211", "Doublet"];

    basic_signals_plot = figure;
    comb_signals_plot = figure;
    counter_comb = 1;
    counter_basic = 1;
    for id = 1:7

        rng default;

        % Assign parameters required for the experiment
        [optimal_input{id}, signal_type{id}, properties] = LoadOptimalInput(id, sample_time);
        amplitudes_cell{id} = properties.amplitudes;

        metric_selector = 1; % Trace of C
        display_flag = 0; % No display

        if isfield(properties, 'tfrac') % For combination with specified tfrac
            [J_cell(id), fitmodel_cell{id}, est_params{id}, est_unc{id}, dataset{id}, excitation_cell{id}] = IdentificationExperiment(optimal_input{id},...
                theta0, signal_type{id}, metric_selector, display_flag, amplitudes_cell{id}, properties.tfrac);
        else
            [J_cell(id), fitmodel_cell{id}, est_params{id}, est_unc{id}, dataset{id}, excitation_cell{id}] = IdentificationExperiment(optimal_input{id},...
                theta0, signal_type{id}, metric_selector, display_flag, amplitudes_cell{id});
        end

        est_err(:, id) = abs(th_true' - est_params{id});


        if id == 1 || id == 2 || id == 4

            figure(basic_signals_plot);
            subplot(3, 1, counter_comb);
            plot(excitation_cell{id}(:, 1), excitation_cell{id}(:, 2), '-', 'LineWidth', 1.05);
            xlabel('Time [s]');
            ylabel('Amplitude [-]');
            title("Signal: " + name_cell(id));
            % Default Options
            grid on;
            axis auto
            ax_gca = gca;
            ax_gca.XAxisLocation = 'bottom';
            ax_gca.YAxisLocation = 'left';
            ax_gca.XMinorTick = 'on';
            ax_gca.YMinorTick = 'on';
            ax_gca.LineWidth = 1.02;
            hold off;

            counter_comb = counter_comb + 1;

        else
            figure(comb_signals_plot);
            subplot(2, 2, counter_basic);
            plot(excitation_cell{id}(:, 1), excitation_cell{id}(:, 2), '-', 'LineWidth', 1.05);
            xlabel('Time [s]');
            ylabel('Amplitude [-]');
            title("Signal: " + name_cell(id));
            % Default Options
            grid on;
            axis auto
            ax_gca = gca;
            ax_gca.XAxisLocation = 'bottom';
            ax_gca.YAxisLocation = 'left';
            ax_gca.XMinorTick = 'on';
            ax_gca.YMinorTick = 'on';
            ax_gca.LineWidth = 1.02;
            hold off;
            counter_basic = counter_basic + 1;

        end

        FIT_perc{id} = fitmodel_cell{id}.Report.Fit.FitPercent;
        disp("Cost for combination: " + num2str(id) + ": " + num2str(J_cell(id)) + "    " + ...
            " FIT: " + num2str(FIT_perc{id}(1)) + "%, " + num2str(FIT_perc{id}(2)) + "%    " + ...
            "Total est. error: " + num2str(sum(est_err(:, id))));

    end

    % Find input with lowest uncertainty
    [minJ, minpos] = min(J_cell);
    sigma = est_unc{minpos};
    estimates = est_params{minpos};

    disp('-------- Best Input results --------')
    fprintf("\n Parameters and relative 3-sigma uncertainty (Gaussian distr. assumption)" + ...
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



%% Identified model validation
%% Simulate reference model
load('input_workspace.mat')
clear ExcitationM

% Generate Validation signal
rng(1);
params.dt = 0.001;
params.f0 = 0;
params.ff = 15;
params.t0 = 0;
params.tf = 150;

[validation_signal, validation_timevec] = GenerateInput(params, 1);
ExcitationM = [validation_timevec, validation_signal];

% Get reference output
sim_obj = SetModel(th_true, ExcitationM);
output = sim(sim_obj);

% Signals Pre-Processing
N_delay = 1;
[Mtot_cell{1}, ax_cell{1}, q_cell{1}, time_grid_cell{1}, pitch_angle_cell{1}] = OutputPreProcess(output, N_delay);

clear output sim_obj

%% Simulate identified model with Task1 input
load('input_workspace.mat')
clear ExcitationM

load('Task1_results.mat')

rng(1);
ExcitationM = [validation_timevec, validation_signal];
sim_object = SetModel(theta1, ExcitationM);
output = sim(sim_object);

J_task1 = sum(est_unc1.^2);

% Signals Pre-Processing
[Mtot_cell{2}, ax_cell{2}, q_cell{2}, time_grid_cell{2}, pitch_angle_cell{2}] = OutputPreProcess(output, 1);

disp('-------- Task 1 Input results --------')
fprintf("\n Parameters and relative 3-sigma uncertainty (Gaussian distr. assumption)" + ...
    "\nXu = %3.3g\t 3sig%%: %3.3g" + ...
    "\nMu = %3.3g\t 3sig%%: %3.3g" + ...
    "\nXq = %3.3g\t 3sig%%: %3.3g" + ...
    "\nMq = %3.3g\t 3sig%%: %3.3g" + ...
    "\nXd = %3.3g\t 3sig%%: %3.3g" + ...
    "\nMd = %3.3g\t 3sig%%: %3.3g", ...
    theta1(1), 100*abs(3*est_unc1(1)./theta1(1)), ...
    theta1(2), 100*abs(3*est_unc1(2)./theta1(2)), ...
    theta1(3), 100*abs(3*est_unc1(3)./theta1(3)), ...
    theta1(4), 100*abs(3*est_unc1(4)./theta1(4)), ...
    theta1(5), 100*abs(3*est_unc1(5)./theta1(5)), ...
    theta1(6), 100*abs(3*est_unc1(6)./theta1(6)));


fprintf('\n\n');
fprintf('-------- Errors with respect to "truth" model --------');

bias_perc1 = 100 * abs((theta1' - th_true)./th_true);

cell_param = {'Xu', 'Mu', 'Xq', 'Mq', 'Xd', 'Md'};

fprintf("\n");
for idp = 1:length(est_unc1)
    fprintf("%s %% estimation error: %4.3g \n",...
        cell_param{idp}, bias_perc1(idp));
end

%% Simulate identified model with optimized input
load('input_workspace.mat')
decimation = 1;
clear ExcitationM

% Hard-coded parameters vector
theta = estimates;

rng(1);
ExcitationM = [validation_timevec, validation_signal];
sim_object = SetModel(theta, ExcitationM);

output = sim(sim_object);

% Signals Pre-Processing
N_delay = 1;

[Mtot_cell{3}, ax_cell{3}, q_cell{3}, time_grid_cell{3}, pitch_angle_cell{3}] = OutputPreProcess(output, 1);


% Evaluate time-domain errors with respect to reference
% Ideally: zero mean signal with random noise due to 
% Either time "error" signal in logarithmic scale, use compare.
for i = 2:3
    output_error{i-1, 1} = q_cell{i} - q_cell{1};
    output_error{i-1, 2} = ax_cell{i} - ax_cell{1};
end

%% Plots
% Frequency Response with bode plots
bode_fig = figure;
w_axis = 0:0.05:(2*pi*100);
bplot = bodeplot(Hmodelstruct(th_true), 'r-',...
    fitmodel_cell{minpos}, 'k-', model_task1, 'b-', w_axis);
showConfidence(bplot, 3);
legend('"Truth" model', 'Identified "optimal" model', 'Identified model task1');
title('Bode diagrams comparison with 3$\sigma$ confidence interval', 'Interpreter', 'Latex')
% Default options applied to all axes handles
axes_handles = findall(bode_fig, 'type', 'axes');
set(axes_handles, 'YMinorTick', 'on', 'XMinorTick', 'on', 'LineWidth', 1.04, ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on');

% Pole-Zero plot
pfig = figure;
pz_plot = iopzplot(Hmodelstruct(th_true), 'r', fitmodel_cell{minpos}, 'k', model_task1, 'b');
showConfidence(pz_plot, 3);

legend('"Truth" model', 'Identified "optimal" model', 'Identified model task1', '', '');
% Default Options
% Default options applied to all axes handles
axes_handles = findall(pfig, 'type', 'axes');
set(axes_handles, 'YMinorTick', 'on', 'XMinorTick', 'on', 'LineWidth', 1.04, ...
    'XMinorGrid', 'on', 'YMinorGrid', 'on');
title('Poles/Zeros comparison with 3$\sigma$ confidence interval', 'Interpreter', 'Latex');

%% Computation of Validation Indexes
% evalFIT = @(qmeas, q) 100 * max([0, 1 - (norm(qmeas - q)).^2 ./ norm(qmeas - mean(qmeas))]);
% evalPEC = @(qmeas, q) 1/sqrt(length(qmeas)) .* (norm(qmeas - q)).^2;
% evalVAF = @(qmeas, q) 100 * max([0, 1 - var(qmeas - q)./var(qmeas)]);
% 
% ModelFIT = [evalFIT(pitch_angle_ref, pitch_angle); evalFIT(q_ref, q); evalFIT(ax_ref, ax)];
% ModelPEC = [evalPEC(pitch_angle_ref, pitch_angle); evalPEC(q_ref, q); evalPEC(ax_ref, ax)];
% ModelVAF = [evalVAF(pitch_angle_ref, pitch_angle); evalVAF(q_ref, q); evalVAF(ax_ref, ax)];

% Display results:
% fprintf("\n\nModelFIT: theta: %3.1g%%, q: %3.1g%%, ax: %3.1g%%\n  ", ModelFIT(1), ModelFIT(2), ModelFIT(3));
% fprintf("ModelPEC: theta: %3.1g, q: %3.1g, ax: %3.1g\n  ", ModelPEC(1), ModelPEC(2), ModelPEC(3));
% fprintf("ModelVAF: theta: %3.1g%%, q: %3.1g%%, ax: %3.1g%%\n  ", ModelVAF(1), ModelVAF(2), ModelVAF(3));



end










