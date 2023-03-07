function [optimal_input, signal_type, params] = LoadOptimalInput(comb, sample_time, optimal_input)
%% PROTOTYPE
% [optimal_input, signal_type, params] = LoadOptimalInput(comb, sample_time, optimal_input)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function storing the best combination found for the parameters defining
% the inputs to the system plant
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% comb: [1] Number identifying the combination: 
%           1) LSW + RBS, 2) LSW + 3211, 3) Sine Sweep only, 4) Linear Sine Sweep + 3211 (80/20)
%           5) RBS, 6) 3211, 7) Doublet
% sample_time: [1] sample time to generate the signal
% optimal_input: [variable number]  optimal combination of parameters as
%                returned by the optimization cycle
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% optimal_input: input signal defined by the parameters
% signal_type: [1] number identifying the signal type as defined in GenerateInput()
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-03-2023    Pietro Califano     Function documented
% -------------------------------------------------------------------------------------------------------------


%% Function code

switch comb

    case 1 % LSW + RBS
        if ~exist('optimal_input', 'var')
            % Achieved J =
            optimal_input = [12.2709802456560, 297.115680758792, 0.144765478881002,	1.30758171976682, 3.10581101055947];
        end

        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        f0 = optimal_input(3);
        ff = optimal_input(4);
        T = optimal_input(5);

        signal_type = [1, 3];


        params.tfrac = 0.5; % Percentage of time for the first signal type
        params.K = K;
        params.t0 = t0;
        params.tf = tfin;
        params.f0 = f0;
        params.ff = ff;
        params.dt = sample_time;
        params.T = T;
        params.amplitudes = [1, 0.5];

    case 2 % LSW + 3211
        if ~exist('optimal_input', 'var')
            % Achieved J =
            optimal_input = [37.6897636715143	259.196074597491	0.196093750000000	1.00008907914162	1.00009547173977];
        end

        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        f0 = optimal_input(3);
        ff = optimal_input(4);
        N = optimal_input(5);

        signal_type = [1, 4];
        params.tfrac = 0.5;  % Percentage of time for the first signal type
        params.K = K;
        params.t0 = t0;
        params.tf = tfin;
        params.f0 = f0;
        params.ff = ff;
        params.dt = sample_time;
        params.N = N;

        params.amplitudes = [1, 0.75];

    case 3 % Sine Sweep only
        % Pameters: 1) K, 2) tf, 3) f0, 4) ff
        if ~exist('optimal_input', 'var')
            % Achieved J = 0.04486
            optimal_input = [31.085075837497065, 2.564445630905003e+02, 0.195833893040528, 1];
        end

        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        f0 = optimal_input(3);
        ff = optimal_input(4);

        signal_type = 1;
        params.K = K;
        params.t0 = t0;
        params.tf = tfin;
        params.f0 = f0;
        params.ff = ff;
        params.dt = sample_time;

        params.amplitudes = 1;

    case 4 % Linear Sine Sweep + 3211 (80/20)
        if ~exist('optimal_input', 'var')
            % Achieved J =
            optimal_input = [20.3173576999760	247.000135224483	0.0471025295387743	0.800061035156250	1];
        end

        % Pameters: 1) K, 2) tf, 3) f0, 4) ff, 5) N
        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        f0 = optimal_input(3);
        ff = optimal_input(4);
        N = optimal_input(5);

        signal_type = [1, 4];
        % Percentage of time for the first signal type
        params.tfrac = 0.8;
        params.K = K;
        params.t0 = t0;
        params.tf = tfin;
        params.f0 = f0;
        params.ff = ff;
        params.dt = sample_time;
        params.N = N;
        params.amplitudes = [1, 0.75];

    case 5 % RBS
        % Pameters: 1) K, 2) tf, 3) T
        if ~exist('optimal_input', 'var')
            % Achieved J =
            optimal_input = [14.3057893512984	295.401715933268	9.70215092871663];
        end

        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        T = optimal_input(3);

        signal_type = 3;
        params.K = K;
        params.t0 = t0;
        params.tf = tfin;
        params.dt = sample_time;
        params.T = T;
        params.amplitudes = 0.5;

    case 6 % 3211
        % Pameters: 1) K, 2) tf, 3) N
        if ~exist('optimal_input', 'var')
            % Achieved J = 0.1079
            optimal_input = [21.751767024881480, 2.842914704650040e+02, 2.368796254243542];
        end
        % NOTE: N is rounded before being used

        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        N = optimal_input(3);

        signal_type = 4;
        params.K = K;
        params.t0 = t0;
        params.tf = tfin;
        params.dt = sample_time;
        params.N = N;
        params.amplitudes = 0.75;

    case 7 % Doublet
        % Pameters: 1) K, 2) tf, 3) N
        if ~exist('optimal_input', 'var')
            % Achieved J = 0.05704
            optimal_input = [14.928329041074418, 2.714867292872423e+02, 2.239365007274643];
        end

        % NOTE: N is rounded before being used

        t0 = 0;
        K = optimal_input(1);
        tfin = optimal_input(2);
        N = optimal_input(3);

        signal_type = 5;
        params.K = K;
        params.t0 = t0;
        params.tf = tfin;
        params.dt = sample_time;
        params.N = N;
        params.amplitudes = 0.75;
end

end