function [fitmodel, est_params, est_unc] = greyest_wrapper(data_to_fit, model_fun, theta0, display_flag)
%% PROTOTYPE
% [fitmodel, est_params, est_unc] = greyest_wrapper(data_to_fit, model_fun, theta0, display_flag)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Wrapper function to call greyest()
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% data_to_fit: [iddata object] dataset passed to greyest() for the identification
% model_fun: function returning the state space model of the system
% theta0: [6x1] initial guess of the parameters vector
% display_flag: [bool] 0: Do not show iterations, 1: Show iterations
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% fitmodel: [idgrey object] fitted model returned by greyest()
% est_params: [6x1] estimated parameters vector
% est_unc: [6x1] 1sigma standard deviation expressing the estimate uncertainty
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-03-2023    Pietro Califano     Function documented
% -------------------------------------------------------------------------------------------------------------


%% Function code
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', 16);

% Guess parameters vector
% theta0 = th_true.*[1.2 1.1 0.8 1.1 0.7 1.1];

% Generate initial guess for greyest function
% odefun = 'LongDyn_ODE';
% Derivative of udot wrt u
Xu = theta0(1);
% Derivative of udot wrt q
Xq = theta0(2);
% Derivative of qdot wrt u
Mu = theta0(3);
% Derivative of qdot wrt q
Mq = theta0(4);
% Derivative of udot wrt delta
Xd = theta0(5);
% Derivative of qdot wrt delta
Md = theta0(6);

% Model parameters cell
parameters = {'Xu', Xu ; 'Xq', Xq; 'Mu', Mu;...
    'Mq', Mq; 'Xd', Xd; 'Md', Md};
% Select continuous-time model
fcn_type = 'c';

% Assembly Frequency Response Data
% data = frd(freqdata, faxis_masked);
% Create grey model for identification
greyobj = idgrey(model_fun, parameters, fcn_type, 'OutputUnit', {'rad/s', 'm/s2'},...
    'InputName', {'M_{tot}'}, 'OutputName', {'q','ax'});

% Set grey identification options
if display_flag == 1
    greyopt = greyestOptions('Display', 'on');
else
    greyopt = greyestOptions('Display', 'off');
end
% Identify model parameters with MATLAB Built-In function greyest()
[fitmodel, ~] = greyest(data_to_fit, greyobj, greyopt);
[est_params, est_unc] = getpvec(fitmodel);

end