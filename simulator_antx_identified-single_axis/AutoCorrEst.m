function [R_xx, signal_zm] = AutoCorrEst(signal, Nsamples)
%% PROTOTYPE
% [R_xx, signal_zm] = AutoCorrEst(signal, Nsamples);
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
x = signal;
R_xx = zeros(Nsamples, 1);

signal_zm = mean(x);
x = x - signal_zm;

for nt = 1:Nsamples-1
    for n = 1:Nsamples - nt
        R_xx(nt) = R_xx(nt) + (x(n)*x(n + nt))/(Nsamples-nt);
    end
end


end