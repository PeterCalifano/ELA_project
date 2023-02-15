function [signal, taxis] = GenerateInput(params, signal_type)
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

        if nargout > 1
            taxis = timevec;
        end


    case 2 % Pseudo Random Binary Sequence
        t0 = params.t0;
        T = params.T;
        dt = T/10;
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

        if nargout > 1
            taxis = timevec;
        end
end