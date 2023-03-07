function [gamma2] = EstimateCoherence(uy, uu, yy)
%% PROTOTYPE
% [gamma2] = EstimateCoherence(uy, uu, yy)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Estimates the coherence function from input, output autospectra and
% cross-spectrum and plots the result.
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% uy: Input/Output cross-spectrum
% uu: Input Autospectrum
% yy: Output Autospectrum 
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% gamma2: estimated coherence of uu and yy
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-03-2023    Pietro Califano     Function documented
% -------------------------------------------------------------------------------------------------------------


%% Function code


Nf = length(uy);
gamma2 = abs(uy).^2 ./ abs(uu.*yy);

% Check that gamma2 did not result to be greater than 1 at any frequency
if sum(gamma2) > Nf
    warning('Coherence at some frequencies resulted greater than one!')
end



end