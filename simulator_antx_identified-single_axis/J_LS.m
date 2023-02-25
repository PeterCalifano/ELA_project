function [J, eH, R] = J_LS(yH, yH_sim, R)
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

[Nf2, Nfcn] = size(yH);
Nf = Nf2/2;
Nfcn2 = 2*Nfcn;

flag_R = 0;
% Default Variance matrix = identity
if nargin < 3
    if (Nf2 * Nfcn2)^2 < 100^2 % Use vectorial operation
        R = eye(Nf2 * Nfcn2);
    else % Use for loop
        R = eye(Nfcn);
    end
else
    flag_R = 1;
end

% Evaluate deviation and cost function at each frequency
J = 0;

% J = 1/2 * sum(e'* (R^-1) * e);

% Compute deviations
eH = reshape(yH - yH_sim, Nf, 2*Nfcn);
% Compute scalar cost
for idf = 1:Nf
    R = eH' * eH;
end
R = R./Nf;
Rinv = R^(-1);

for id = 1:Nf
    J = J + 0.5 * eH(id, :) * Rinv * eH(id, :)';
end

end