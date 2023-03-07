function [signal, timevec] = CombineInput(params, signal_type, amplitude)
%% PROTOTYPE
% [signal, timevec] = CombineInput(params, signal_type, amplitude)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Wrapper of GenerateInput()
% Generates input signals and combines them (up to 2 different types) such
% that the total signal has duration tf, of which the first covers tfrac % 
% of the time. For Linear Sine Sweep, f0 and ff determines the bandwith. T
% determines the minimum time for which a level is maintained in the RBS, N
% is the number of repetitions for 3211 and Doublet.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% params: [struct] with variable fields: t0, tf, f0, ff, dt, T, N, tfrac. 
% signal_type: scalar or [2x1] vector specifying the type of signal
%              1) Linear Sine Sweep, 2) Logarithmic Sine Sweep
%              3) Random Binary Sequence, 4) 3211 sequence, 5) Doublet
% amplitude: scalar or [2x1] specifying the scaling factor of the
%            amplitude. (default = 1)
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% signal: [N_instants x 1] generated signal, with length given as
% N_instants = ceil(tf - t0)/dt 
% timevec: [N_instants x 1] timegrid on which the signal is generated
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-03-2023    Pietro Califano     Function documented
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% GenerateInput() function
% -------------------------------------------------------------------------------------------------------------

%% Function code
% Check for fraction in params
if isfield(params, 'tfrac') && length(signal_type) > 1
    tfrac = params.tfrac;
elseif length(signal_type) > 1
    tfrac = 0.5;
else
    tfrac = 1;
end

if nargin < 3
   amplitude = ones(length(signal_type), 1);
end

tf_mid = tfrac * params.tf;
counter = 1;
N_instants = ceil(params.tf - params.t0/params.dt);

timevec = nan(N_instants, 1);
signal = nan(N_instants, 1);

% Allocation indexes
idfirst = 1;
idlast = 0;

for input_type = signal_type

    params_temp = params;
    % Assign tf_mid for first signal in place of tf
    if counter == 1 && length(signal_type) > 1
        params_temp.tf = tf_mid;
    elseif counter == 2 
        params_temp.t0 = tf_mid + params.dt;
    end

    [signal_temp, time_temp] = GenerateInput(params_temp, input_type);
   
    idlast = idlast + length(signal_temp);

    timevec(idfirst:idlast, 1) = time_temp;
    signal(idfirst:idlast, 1) = amplitude(counter).*signal_temp;

    idfirst = idlast + 1;

    counter = counter + 1;
end

signal = rmmissing(signal);
timevec = rmmissing(timevec);

end