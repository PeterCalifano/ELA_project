% Optimization of the experiment for better model identification
theta0 = th_true.*[1.2 1.1 0.8 1.1 0.7 1.1];
% x0 = [8, 295, 0.02, 50];
x0 = [8,295,20];
metric_selector = 1;


J = IdentificationExperiment(x0, theta0, metric_selector);
disp(J);


fmincon_opts = optimoptions("fmincon", "Algorithm", 'interior-point',...
       'FunctionTolerance', 1e-5, 'MaxIterations', 10,...
       'UseParallel', false, 'Display', 'iter', 'OptimalityTolerance', 1e-6);

ga_opts = optimoptions("ga", "Display", 'iter', 'CrossoverFraction', 0.7, ...
    'FunctionTolerance', 1e-4, 'MaxTime', 10*60, 'PopulationSize', 100, 'UseParallel', false);

% sine sweep
% LB = [5, 1*60, 0.001, 1]; 
% UB = [40, 5*60, 0.1, 150]; 

% 3211
LB = [5, 1*60, 5]; 
UB = [40, 5*60, 50]; 

% [optimal_input, ] = fmincon(@(x) IdentificationExperiment(x, theta0,...
%     metric_selector), x0, [], [], [], [], LB, UB, [], fmincon_opts);

% 3211 [8.173686814228478,2.521058987690398e+02,20.019986773376180]

% [optimal_input, ] = ga(@(x) IdentificationExperiment(x, theta0,...
%     metric_selector), length(LB), ga_opts);


x = [8.173686814228478,2.521058987690398e+02,20.019986773376180];

% Unpack parameters
t0 = 0;
K = x(1);
tfin = x(2);
% T = x(2);
% f0 = x(3);
% ff = x(4);

 N = x(3);

signal_type = 4;

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
%% Local Functions
function J = IdentificationExperiment(x, theta0, metric_selector)

% Load input (reference) model
load('input_workspace.mat')
clear ExcitationM;

% Unpack parameters
t0 = 0;
K = x(1);
tf = x(2);
% T = x(2);
% f0 = x(3);
% ff = x(4);

 N = x(3);

signal_type = 4;

% Generate Input
% params.t0 = t0;
% params.tf = tf;
% params.T = T;

% params.t0 = t0;
params.tf = tf;
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

[~, q_zm, ~] = AutoCorrEst(q, Nsamples);
[~, ax_zm, ~] = AutoCorrEst(ax, Nsamples);
[~, delta_zm, ~] = AutoCorrEst(Mtot, Nsamples);

N = length(delta_zm);

Window = hanning(floor(N/K));

[cohr1, f1] = mscohere(delta_zm, q_zm, Window, [], [], 1/sample_time); % f1 is in Hz
[cohr2, ~] = mscohere(delta_zm, ax_zm, Window, [], [], 1/sample_time);

f_window = f1(cohr1 > 0.6 & cohr2 > 0.6);
% Define data structure for greyest
data_time = iddata([q_zm, ax_zm], delta_zm, sample_time);
% Transform to frequency domain
data = fft(data_time); % frequency in rad/s

q_f = data.OutputData(:, 1);
ax_f = data.OutputData(:, 2);
delta_f = data.InputData;
faxis = data.Frequency ./ (2*pi);

f_lb = f_window(1);
f_ub = f_window(end);

id_f = faxis >= f_lb & faxis <= f_ub;
q_data = q_f(id_f);
ax_data = ax_f(id_f);
delta_data = delta_f(id_f);

data_to_fit = iddata([q_data, ax_data], delta_data, sample_time,...
    'Frequency', faxis(id_f), 'FrequencyUnit', 'Hz');

% Generate initial guess for greyest function
model_fun = 'LongDyn_ODE';
[~, est_params, est_unc] = greyest_wrapper(data_to_fit, model_fun, theta0);

% Compute Objective function Cost

switch metric_selector
    case 1
        J = sum(est_unc);
    case 2
        J = max(est_unc);
    case 3
        J = det(diag(est_unc));
    case 4
        J = sum(est_unc./est_params);
end


end

