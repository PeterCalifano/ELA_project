function [J, fitmodel, est_params, est_unc, dataset, excitation] = IdentificationExperiment(x, theta0, signal_type, metric_selector, display_flag, amplitude, tfrac)

% Load input (reference) model
load('input_workspace.mat')
clear ExcitationM;

if ~exist('amplitude', 'var')
   amplitude = ones(length(signal_type), 1);
end

% Unpack parameters
t0 = 0;
K = x(1);
tf = x(2);

if exist('tfrac', 'var')
    params.tfrac = tfrac;
end

params.t0 = t0;
params.tf = tf;
params.dt = sample_time;

if length(signal_type) > 1 || signal_type <= 2
    f0 = x(3);
    ff = x(4);
    params.ff = ff;
    params.f0 = f0;
end

% Generate Input
if length(signal_type) > 1
    switch signal_type(end)
        case 3
            T = x(5);
            params.T = T;
        case 4
            N = x(5);
            params.N = N;
        case 5
            N = x(5);
            params.N = N;
    end

elseif signal_type == 3
    T = x(3);
    params.T = T;
elseif signal_type == 4 || signal_type == 5
    N = x(3);
    params.N = N;
end


try

    [signal, timevec] = CombineInput(params, signal_type, amplitude);
    ExcitationM = [timevec, signal];

    % Time grid
    t = ExcitationM(:, 1);
    simulation_time = t(end) - t(1);

    % Simulate system
    sim_object = Simulink.SimulationInput(model_name);
    sim_object = setVariable(sim_object, 'ExcitationM', ExcitationM);
    sim_object = sim_object.setModelParameter('StopTime', num2str(simulation_time));
    sim_object = sim_object.setModelParameter('SimulationMode', 'accelerator');

    output = sim(sim_object);

    % Signals Pre-Processing
    N_delay = 1;

    % Extract useful input/output samples
    Mtot = output.Mtot;
    ax = output.ax;
    q = output.q;

    % dt = 1/250; % 250 Hz, defined in parameters_controller
    % Consider delay of the output (4 samples)
    Mtot = Mtot(1:(end-N_delay));
    ax = ax((1+N_delay):end);
    q = q((1+N_delay):end);

    %% Model identification process
    % Remove Mean from signals
    q_zm = q - mean(q);
    ax_zm = ax - mean(ax);
    delta_zm = Mtot - mean(Mtot);

    N_T = length(delta_zm);

    Window = hanning(floor(N_T/K));

    [cohr1, f1] = mscohere(delta_zm, q_zm, Window, [], [], 1/sample_time); % f1 is in Hz
    [cohr2, ~] = mscohere(delta_zm, ax_zm, Window, [], [], 1/sample_time);

    f_window = f1(cohr1 >= 0.6 & cohr2 >= 0.6);
    % Define data structure for greyest
    data_time = iddata([q_zm, ax_zm], delta_zm, sample_time);
    % Transform to frequency domain
    data = fft(data_time); 

    q_f = data.OutputData(:, 1);
    ax_f = data.OutputData(:, 2);
    delta_f = data.InputData;
    % Frequency in Hz
    faxis = data.Frequency ./ (2*pi);

    if length(f_window) < 20 % Break cycle if window of available frequency is too narrow
        warning('Frequency range with less than 20 samples');
        J = 1e10;
        return;
    end

    f_lb = f_window(1);
    f_window = f_window(f_window <= 15);
    f_ub = f_window(end);

    id_f = faxis >= f_lb & faxis <= f_ub;
    q_data = q_f(id_f);
    ax_data = ax_f(id_f);
    delta_data = delta_f(id_f);

    data_to_fit = iddata([q_data, ax_data], delta_data, sample_time,...
        'Frequency', faxis(id_f), 'FrequencyUnit', 'Hz');

    % Generate initial guess for greyest function
    model_fun = 'LongDyn_ODE';
    [fitmodel, est_params, est_unc] = greyest_wrapper(data_to_fit, model_fun, theta0, display_flag);

    % Create output  if needed
    if nargout > 4
        dataset.data_time = data_time;
        dataset.output = output;
        dataset.data_to_fit = data_to_fit;
        if nargout > 5
            time_grid = output.time_grid;
            Excit_signal = output.Excit_signal(1: (end-N_delay) );
            excitation = [time_grid((1+N_delay):end), Excit_signal];
        end
    end

    % Compute Objective function Cost
    switch metric_selector
        case 1
            J = sum(est_unc.^2);
        case 2
            J = max(est_unc.^2);
        case 3
            J = det(diag(est_unc.^2));
        case 4
            J = sum((est_unc./est_params).^2);   
    end

catch
    J = 1e10;
end


end
