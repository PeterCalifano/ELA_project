% Optimization of the experiment for better model identification
theta0 = th_true.*[1.2 1.1 0.8 1.1 0.7 1.1];
% theta0 = -rand(1, 6).*th_true;

comb = 2;

switch comb
    case 1 % Linear Sine Sweep + RBS
        % Pameters: 1) K, 2) tf, 3) f0, 4) ff, 5) T

        % Initial guess for optimization
        x0 = [8, 295, 0.05, 20, 1];

        LB = [5, 1*60, 0.001, 1, 0.05];
        UB = [40, 5*60, 0.2, 150, 10];

        signal_type = [1, 3];

    case 2 % Linear Sine Sweep + 3211
        % Pameters needed
        % Pameters: 1) K, 2) tf, 3) f0, 4) ff, 5) N

        % Initial guess for optimization
        x0 = [8, 295, 0.05, 20, 6];

        LB = [5, 1*60, 0.001, 1, 1];
        UB = [40, 5*60, 0.2, 150, 20];
        signal_type = [1, 4];

    case 3 % Sine Sweep 
        % Pameters needed
        % Pameters: 1) K, 2) tf, 3) f0, 4) ff

        % Initial guess for optimization
        x0 = [8, 295, 0.05, 20];

        LB = [5, 1*60, 0.001, 1];
        UB = [40, 5*60, 0.2, 150];

        signal_type = 1;

    case 4 % Logarithmic Sine Sweep
        % Pameters needed
        % Pameters: 1) K, 2) tf, 3) f0, 4) beta

        % Initial guess for optimization
        x0 = [8, 295, 0.05];

        LB = [5, 1*60, 0.001];
        UB = [40, 5*60, 0.2];

        signal_type = 2;

    case 5 % RBS

        % Pameters needed
        % Pameters: 1) K, 2) tf, 3) T

        % Initial guess for optimization
        x0 = [8, 295, 1];

        LB = [5, 1*60, 0.01];
        UB = [40, 5*60, 10];

        signal_type = 3;

    case 6 % 3211

        % Pameters needed
        % Pameters: 1) K, 2) tf, 3) N

        % Initial guess for optimization
        x0 = [8, 295, 5];

        LB = [5, 1*60, 2];
        UB = [40, 5*60, 30];

        signal_type = 4;

    case 7 % Doublet
        % Pameters needed
        % Pameters: 1) K, 2) tf, 3) N

        % Initial guess for optimization
        x0 = [8, 295, 5];

        LB = [5, 1*60, 2];
        UB = [40, 5*60, 50];

        signal_type = 5;

end

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
    'MaxStallGenerations', 10,'UseParallel', false);

% Run optimization
% [optimal_input, ] = fmincon(@(x) IdentificationExperiment(x, theta0,...
%     metric_selector, 0), x0, [], [], [], [], LB, UB, [], fmincon_opts);

% 3211 [8.173686814228478,2.521058987690398e+02,20.019986773376180]

[optimal_input, ] = ga(@(x) IdentificationExperiment(x, theta0, signal_type,...
    metric_selector, 0), length(LB), [], [], [], [], LB, UB, [], [], ga_opts);


[optimal_input2, ] = fmincon(@(x) IdentificationExperiment(x, theta0, signal_type,...
    metric_selector, 0), optimal_input, [], [], [], [], LB, UB, [], fmincon_opts);


save('comb2_gafmincon.mat');

return


%%

% x = [8.173686814228478,2.521058987690398e+02,20.019986773376180];
% comb = 3
% Unpack parameters
switch comb

    case 1 % LSW + RBS
        t0 = 0;
        K = optimal_input2(1);
        tfin = optimal_input2(2);
        f0 = optimal_input2(3);
        ff = optimal_input2(4);
        T = optimal_input2(5);

        params.T = T;

    case 2 % LSW + 3211

        t0 = 0;
        K = optimal_input2(1);
        tfin = optimal_input2(2);
        f0 = optimal_input2(3);
        ff = optimal_input2(4);
        N = optimal_input2(5);

        params.N = N;

end

% T = x(2);

% signal_type = 4;

% Generate Input
% params.t0 = t0;
% params.tf = tf;
% params.T = T;


% Create params to test and retrieve theta vector
params.ff = ff;
params.f0 = f0;

params.t0 = t0;
params.tf = tfin;
params.dt = sample_time;


[signal, timevec] = GenerateInput(params, signal_type);
ExcitationM = [timevec, signal];

% Time grid
t = ExcitationM(:, 1);
simulation_time = t(end) - t(1);

% Simulate system
sim_object = Simulink.SimulationInput(model_name);
sim_object = setVariable(sim_object, 'ExcitationM', ExcitationM);
sim_object = sim_object.setModelParameter('StopTime', num2str(simulation_time));

output = sim(sim_object);

% Signals Pre-Processing
N_delay = 1;

Excit_signal = output.Excit_signal;

% Extract useful input/output samples
Mtot = output.Mtot;
time_grid = output.time_grid;
ax = output.ax;
q = output.q;

% dt = 1/250; % 250 Hz, defined in parameters_controller
time_grid = time_grid((1+N_delay):end);
% Consider delay of the output (4 samples)
Mtot = Mtot(1:(end-N_delay));
ax = ax((1+N_delay):end);
q = q((1+N_delay):end);

% Model identification process
Nsamples = length(time_grid);

IdentifyModel;

fprintf("\nIdentified model parameters and relative 3-sigma uncertainty (Gaussian distr. assumption):" + ...
    "\nXu = %3.3g\t 3sig\%: %3.3g" + ...
    "\nMu = %3.3g\t 3sig\%: %3.3g" + ...
    "\nXq = %3.3g\t 3sig\%: %3.3g" + ...
    "\nMq = %3.3g\t 3sig\%: %3.3g" + ...
    "\nXd = %3.3g\t 3sig\%: %3.3g" + ...
    "\nMd = %3.3g\t 3sig\%: %3.3g", ...
    theta(1), 100*abs(3*est_unc(1)./theta(1)), ...
    theta(2), 100*abs(3*est_unc(2)./theta(2)), ...
    theta(3), 100*abs(3*est_unc(3)./theta(3)), ...
    theta(4), 100*abs(3*est_unc(4)./theta(4)), ...
    theta(5), 100*abs(3*est_unc(5)./theta(5)), ...
    theta(6), 100*abs(3*est_unc(6)./theta(6)));

%% Identified model validation
% Generate input signal for validation

params.ff = 3;
params.f0 = 0.001;

params.t0 = 0;
params.tf = 100;
params.dt = 4e-3;

[validation_signal, validation_timevec] = GenerateInput(params, 1);
validation_signal = validation_signal;
% [validation_signal, validation_timevec] = CombineInput();

%% Simulate reference model
load('input_workspace.mat')
decimation = 1;
clear ExcitationM

rng(1);

ExcitationM = [validation_timevec, validation_signal];

% Time grid
t = ExcitationM(:, 1);
simulation_time = t(end) - t(1);

% Simulate system
sim_object = Simulink.SimulationInput(model_name);
sim_object = setVariable(sim_object, 'ExcitationM', ExcitationM);
sim_object = sim_object.setModelParameter('StopTime', num2str(simulation_time));

output = sim(sim_object);

% Signals Pre-Processing
N_delay = 1;

[Mtot, ax, q, timegrid] = OutputPreProcess(output, 1);

Mtot_cell{1} = Mtot;
ax_cell{1} = ax;
q_cell{1} = q;
time_grid_cell{1} = timegrid;

%% Simulate identified model with Task1 input
load('input_workspace.mat')
clear ExcitationM

% Hard-coded parameters vector
theta = [-2.464780033101573e-01;
     1.322467636559020e-01;
    -5.331283778077723e+00;
    -3.017675910479900e+00;
    -6.624738689764718e+00;
     4.570041745162848e+02];

rng(1);

ExcitationM = [validation_timevec, validation_signal];

sim_object = SetModel(theta, ExcitationM);

output = sim(sim_object);

% Signals Pre-Processing
N_delay = 1;

[Mtot_cell{2}, ax_cell{2}, q_cell{2}, time_grid_cell{2}] = OutputPreProcess(output, 1);

%% Simulate identified model with optimized input
load('input_workspace.mat')
clear ExcitationM

% Hard-coded parameters vector
theta = th_true.*0.001;

rng(1);

ExcitationM = [validation_timevec, validation_signal];

sim_object = SetModel(theta, ExcitationM);

output = sim(sim_object);

% Signals Pre-Processing
N_delay = 1;

[Mtot_cell{3}, ax_cell{3}, q_cell{3}, time_grid_cell{3}] = OutputPreProcess(output, 1);


% Evaluate time-domain errors with respect to reference
% Ideally: zero mean signal with random noise due to 
% Either time "error" signal in logarithmic scale, use compare.
for i = 2:3
    output_error{i-1, 1} = q_cell{i} - q_cell{1};
    output_error{i-1, 2} = ax_cell{i} - ax_cell{1};
end

%% Plots
% Pitch rate
figure;
subplot(2, 1, 1)
hold on;
plot(time_grid_cell{1}, q_cell{1}, '-', 'Linewidth', 1.02, 'DisplayName', 'Reference');
plot(time_grid_cell{2}, q_cell{2}, '-', 'Linewidth', 1.02, 'DisplayName', 'Task1');
plot(time_grid_cell{3}, q_cell{3}, '-', 'Linewidth', 1.02, 'DisplayName', 'Task2');
title('Pitch rate output - Time domain');
grid minor
axis auto;
legend();

subplot(2, 1, 2)
hold on;
plot(time_grid_cell{2}, output_error{1, 1}, '-', 'Linewidth', 1.02, 'DisplayName', 'Task1');
plot(time_grid_cell{3}, output_error{2, 1}, '-',  'Linewidth', 1.02, 'DisplayName', 'Task2');
title('Pitch rate errors - Time domain');
grid minor
axis auto;
legend();

% Acceleration
figure;
subplot(2, 1, 1)
hold on;
plot(time_grid_cell{1}, ax_cell{1}, '-', 'Linewidth', 1.02, 'DisplayName', 'Reference');
plot(time_grid_cell{2}, ax_cell{2}, '-', 'Linewidth', 1.02, 'DisplayName', 'Task1');
plot(time_grid_cell{3}, ax_cell{3}, '-', 'Linewidth', 1.02, 'DisplayName', 'Task2');
title('Acceleration output - Time domain');
grid minor
axis auto;
legend();

subplot(2, 1, 2)
hold on;
plot(time_grid_cell{2}, output_error{1, 2}, '-', 'Linewidth', 1.02, 'DisplayName', 'Task1');
plot(time_grid_cell{2}, output_error{2, 2}, '-', 'Linewidth', 1.02, 'DisplayName', 'Task2');
title('Acceleration errors - Time domain');
grid minor
axis auto;
legend();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










