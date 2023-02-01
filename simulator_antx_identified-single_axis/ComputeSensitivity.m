function dTFdth = ComputeSensitivity(TFfreqresp, freq, params)
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
% Determine size of the parameters vector
Nparam = length(params);
% Static allocation of output
dTFdth = zeros(4*length(freq), Nparam);

% For each parameter theta(idp), compute gradient and store it
for idp = 1:Nparam
    % Reference: MSAS lecture 03
    % Compute perturbation using machine precision value
    dth = max(sqrt(eps), sqrt(eps)*abs(params(idp)));
    % Expand dth to vector
    th_perturbation = zeros(1, length(params));
    th_perturbation(idp) = dth;
    % Determine pertubed theta
    theta = params + th_perturbation;
%     gfreq = squeeze(freqresp(Hmodelstruct(theta), freq));

%     g1re = real(gfreq(1, :));
%     g1im = imag(gfreq(1, :));
% 
%     g2re = real(gfreq(2, :));
%     g2im = imag(gfreq(2, :));

    % "Measurements" vector
%     yH_pert = [g1re, g1im, g2re, g2im]';

    yH_pert = evalFreqR(Hmodelstruct(theta), freq);
    dTFdth(:, idp) = reshape((yH_pert - TFfreqresp)./dth, 4*length(freq), 1);
end

end