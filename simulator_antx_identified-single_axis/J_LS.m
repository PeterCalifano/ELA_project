function [J, eH] = J_LS(yH, yH_sim, R)
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
if 0 %(Nf2 * Nfcn2)^2 < 100^2 && flag_R == 0 % Use vectorial operation
%     warning('TO VALIDATE')
    % Compute deviations
    % 1st col: re(G1), 2nd col: im(G1) to Nth-1 col: re(GN), Nth col: im(GN)
    eH = reshape(yH - yH_sim, Nfcn2*Nf, 1);
    % Compute scalar cost
    J = 0.5 * eH' * R^-1 * eH;

else % Use for loop
    % Compute deviations
    eH = reshape(yH - yH_sim, 2*Nf, Nfcn);
    % Compute scalar cost
%     warning('TO VALIDATE')
    for id = 1:Nf2
        J = J + 0.5 * eH(id, :) * R^-1 * eH(id, :)';
    end

end


end