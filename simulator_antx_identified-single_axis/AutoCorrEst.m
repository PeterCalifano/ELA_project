function [R_xx, signal_zm, signal_mean] = AutoCorrEst(signal, Nsamples)
%% PROTOTYPE
% [R_xx, signal_zm, signal_mean] = AutoCorrEst(signal, Nsamples)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Function estimating the Autocorrelation function of the input signal with
% number of samples given by Nsamples.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% signal [Nsamples x 1]: signals of which to compute the autocorrelation
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% R_xx: [Nsamples x 1]: estimated autocorrelation function
% signal_zm: [Nsamples x 1]: zero mean signal
% signal_mean: [1]: estimated mean of the signal
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-03-2023    Pietro Califano     First version
% -------------------------------------------------------------------------------------------------------------


%% Function code
x = signal;
R_xx = zeros(Nsamples, 1);

% Evaluate mean
signal_mean = mean(x);
% Remove mean
x = x - signal_mean;
% Assign zero mean value as output
signal_zm = x;

% Estimate R_xx
for nt = 1:Nsamples-1
    for n = 1:Nsamples - nt
        R_xx(nt) = R_xx(nt) + (x(n)*x(n + nt))/(Nsamples-nt);
    end
end


end

