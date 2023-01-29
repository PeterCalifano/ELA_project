function [R_xy, R_yx] = CrossCorrEst(signalzm1, signalzm2, Nsamples)
%% PROTOTYPE
% [R_xy, R_yx] = CrossCorrEst(signalzm1, signalzm2, Nsamples)
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
x = signalzm1;
y = signalzm2;

R_xy = zeros(Nsamples, 1);
R_yx = zeros(Nsamples, 1);

% Estimate R_xy
for nt = 1:Nsamples-1
    for n = 1:Nsamples - nt
        R_xy(nt) = R_xy(nt) + (x(n)*y(n + nt))/(Nsamples-nt);
        R_yx(nt) = R_yx(nt) + (y(n)*x(n + nt))/(Nsamples-nt);
    end
end

end