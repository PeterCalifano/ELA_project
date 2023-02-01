function [gamma2, plot_obj] = EstimateCoherence(Gxy, Gxx, Gyy, f_axis)
%% PROTOTYPE
% gamma2 = EstimateCoherence(Gxy, Gxx, Gyy, f_axis)
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
DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

Nf = length(Gxy);

gamma2 = abs(Gxy).^2 ./ (Gxx.*Gyy);

% Check that gamma2 did not result to be greater than 1 at any frequency
if sum(gamma2) > Nf
    warning('Coherence at some frequencies resulted greater than one!')
end

figure;
plot_obj = semilogx(f_axis, gamma2, '.', 'Color', '#aa5533', 'MarkerSize', 4);

xlabel("Frequency $\omega [rad/s]$", 'Interpreter', 'latex');
ylabel("$\gamma_{xy}^2$ [-]", 'Interpreter', 'latex');
title("Estimation of $\gamma_{xy}^2(\omega)$", 'Interpreter', 'latex');
grid minor;
axis auto;


end