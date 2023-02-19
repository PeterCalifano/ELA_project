function [signal, timevec, f_signal, S] = GenerateInput(params, signal_type)
%% PROTOTYPE
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% What the function does
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% in1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% out1 [dim] description
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% Date, User, brief summary of the modification
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades


%% Function code

switch signal_type

    case 1 % Linear Sine Sweep
        t0 = params.t0;
        tf = params.tf;
        ff = params.ff;
        dt = params.dt;
        f0 = params.f0; % [Hz]
        timevec = t0:dt:tf;

        f_signal = f0 + (ff-f0) * (timevec - t0)./(tf-t0);
        signal = sin(2*pi*f_signal.*timevec);


    case 2 % Logarithmic (exp) Sine Sweep

        t0 = params.t0;
        tf = params.tf;
        beta = params.beta;
        dt = params.dt;
        f0 = params.f0; % [Hz]

        timevec = t0:dt:tf;

        f_signal = f0 * exp(beta * timevec);
        signal = sin(2*pi*f_signal.*timevec);

    case 3 % Pseudo Random Binary Sequence
        t0 = params.t0;
        T = params.T;
        dt = T/100;
        tf = params.tf;
        timevec = t0:dt:tf;

        N = ceil(tf/T) + 1;
        RBS = (-1).^(round(rand(N, 1)));

        id = 0;
        signal = zeros(N, 1);

        for i = 1:length(timevec)
            if mod(timevec(i), T) == 0
                id = id + 1;
                signal(i) = RBS(id);        
            else
                signal(i) = RBS(id);
            end

        end

        if nargout > 2
            n = 2^nextpow2(length(timevec));    
            S = fft(signal, n);
            faxis = (1/dt)*(0:(n/2))/n;
            S = S((n/2):end);       
            f_signal = faxis;
        
        end

    case 4 % 3211 sequence
%         t0 = params.t0;
        N = params.N;
        dt = 1/500;
        tf = params.tf;

        TU = tf./(8*N);
        tunit_length = floor(TU./dt); % N° of time instants in 1 time unit
        
        % Build signal_unit
        timeaxis = linspace(0, 8*TU, tunit_length);

        % No input up to 1 TU
        signal_unit(1:tunit_length) = 0; %zeros(tunit_length, 1);
        % HIGH for 3 TU
        signal_unit(tunit_length+1 : 4*tunit_length) = 1;
        % LOW for 2 TU
        signal_unit(4*tunit_length+1 : 6*tunit_length) = -1;
        % HIGH for 1 TU
        signal_unit(6*tunit_length+1 : 7*tunit_length) = 1;
        % LOW for 1 TU
        signal_unit(7*tunit_length+1 : 8*tunit_length) = -1;

        signal_unit = signal_unit';

        % Concatenate N signals
        signal = repmat(signal_unit, N, 1);
        zero_unit = zeros(tunit_length, 1);
        signal = [signal; zero_unit];

        if length(signal) ~= (8*N+1) * tunit_length 
            warning('3211 Signal generation may be wrong');
        end

        if nargout > 1
            timevec = linspace(0, (8*N+1)*TU, length(signal));
        end

    case 5 % Doublet

        N = params.N;
        dt = 1/500;
        tf = params.tf;

        TU = tf./(3*N);
        tunit_length = floor(TU./dt); % N° of time instants in 1 time unit
        
        % Build signal_unit
        timeaxis = linspace(0, 3*TU, tunit_length);

        % No input up to 1 TU
        signal_unit(1:tunit_length) = 0; %zeros(tunit_length, 1);
        % LOW for 1 TU
        signal_unit(tunit_length+1 : 2*tunit_length) = -1;
        % HIGH for 1 TU
        signal_unit(2*tunit_length+1 : 3*tunit_length) = 1;

        signal_unit = signal_unit';

        % Concatenate N signals
        signal = repmat(signal_unit, N, 1);
        zero_unit = zeros(tunit_length, 1);
        signal = [signal; zero_unit];

        if length(signal) ~= (3*N+1) * tunit_length 
            warning('Doublet Signal generation may be wrong');
        end

        if nargout > 1
            timevec = linspace(0, (3*N+1)*TU, length(signal));
        end


end

if isrow(timevec) 
    timevec = timevec';
end

if isrow(signal) 
    signal = signal';
end


end