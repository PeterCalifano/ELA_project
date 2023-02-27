function [signal, timevec] = CombineInput(params, signal_type, amplitude)

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