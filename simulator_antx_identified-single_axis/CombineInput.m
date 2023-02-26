function [signal, timevec] = CombineInput(params, signal_type)

% Check for fraction in params
if isfield(params, 'tfrac')
    tfrac = params.tfrac;
else
    tfrac = 0.5;
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
    signal(idfirst:idlast, 1) = signal_temp;

    idfirst = idlast + 1;

    counter = counter + 1;
end

signal = rmmissing(signal);
timevec = rmmissing(timevec);

end