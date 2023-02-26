function [gamma2, plot_obj] = EstimateCoherence(uy, uu, yy)
%% PROTOTYPE
% [gamma2, plot_obj] = EstimateCoherence(uy, uu, yy, f_axis)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Estimates the coherence function from input, output autospectra and
% cross-spectrum and plots the result.
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


Nf = length(uy);
gamma2 = abs(uy).^2 ./ abs(uu.*yy);

% Check that gamma2 did not result to be greater than 1 at any frequency
if sum(gamma2) > Nf
    warning('Coherence at some frequencies resulted greater than one!')
end



end