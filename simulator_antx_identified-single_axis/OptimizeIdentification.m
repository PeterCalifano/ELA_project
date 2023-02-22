% Optimization of the experiment for better model identification
theta0 = th_true.*[1.2 1.1 0.8 1.1 0.7 1.1];
% theta0 = -rand(1, 6).*th_true;

comb = 1;

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
       'FunctionTolerance', 1e-5, 'MaxIterations', 10,...
       'UseParallel', false, 'Display', 'iter', 'OptimalityTolerance', 1e-6);

ga_opts = optimoptions("ga", "Display", 'iter', 'CrossoverFraction', 0.7, ...
    'FunctionTolerance', 1e-4, 'MaxTime', 5*3600, 'PopulationSize', 100, 'UseParallel', false);

% Run optimization
% [optimal_input, ] = fmincon(@(x) IdentificationExperiment(x, theta0,...
%     metric_selector, 0), x0, [], [], [], [], LB, UB, [], fmincon_opts);

% 3211 [8.173686814228478,2.521058987690398e+02,20.019986773376180]

[optimal_input, ] = ga(@(x) IdentificationExperiment(x, theta0, signal_type,...
    metric_selector, 0), length(LB), [], [], [], [], LB, UB, [], [], ga_opts);

[optimal_input2, ] = fmincon(@(x) IdentificationExperiment(x, theta0, signal_type,...
    metric_selector, 0), optimal_input, [], [], [], [], LB, UB, [], fmincon_opts);

save('comb1_gafmincon.mat');
return


%%

x = [8.173686814228478,2.521058987690398e+02,20.019986773376180];
comb = 3
% Unpack parameters
switch comb

    case 1 % LSW + RBS
        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        f0 = optimal_input(3);
        ff = optimal_input(4);
        T = optimal_input(5);

    case 2 % LSW + 3211

        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        f0 = optimal_input(3);
        ff = optimal_input(4);
        N = optimal_input(5);

end

% T = x(2);

% signal_type = 4;

% Generate Input
% params.t0 = t0;
% params.tf = tf;
% params.T = T;

% params.t0 = t0;
params.tf = tfin;
% params.dt = sample_time;
% params.ff = ff;
% params.f0 = f0;
params.N = N;

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
% RemoveZeroInputMask = Excit_signal ~= 0;

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



% Prototype objective function:
% x: [Nx1] decision variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%