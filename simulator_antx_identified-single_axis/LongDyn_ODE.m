function [A, B, C, D] = LongDyn_ODE(Xu, Xq, Mu, Mq, Xd, Md, Ts)
%% PROTOTYPE
% [A, B, C, D] = LongDyn_ODE(Xu, Xq, Mu, Mq, Xd, Md, Ts)
% -------------------------------------------------------------------------------------------------------------
%% DESCRIPTION
% Builds the state space model of the Longitudinal dynamics
% -------------------------------------------------------------------------------------------------------------
%% INPUT
% Xu: [1] Derivative of udot wrt u
% Xq: [1] Derivative of udot wrt q
% Mu: [1] Derivative of qdot wrt u
% Mq: [1] Derivative of qdot wrt q
% Xd: [1] Derivative of udot wrt delta
% Md: [1] Derivative of qdot wrt delta
% Ts: [1] Sample time
% -------------------------------------------------------------------------------------------------------------
%% OUTPUT
% A: [3x3] Dynamical matrix
% B: [3x1] Input Matrix
% C: [2x3] Output Matrix
% D: [2x1] Feed-through Matrix
% -------------------------------------------------------------------------------------------------------------
%% CHANGELOG
% 07-03-2023    Pietro Califano   Function documented
% -------------------------------------------------------------------------------------------------------------
%% DEPENDENCIES
% -------------------------------------------------------------------------------------------------------------
%% Future upgrades

%% Function code
% Model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)

% Derivative of udot wrt u
% Xu = th(1);
% Derivative of udot wrt q
% Xq = th(2);
% Derivative of qdot wrt u
% Mu = th(3);
% Derivative of qdot wrt q
% Mq = th(4);
% Derivative of udot wrt delta
% Xd = th(5);
% Derivative of qdot wrt delta
% Md = th(6);

% th_true = [Xu, Xq, Mu, Mq, Xd, Md];

%% State space model
A = [Xu, Xq, -9.81; 
    Mu, Mq, 0; 
    0, 1, 0];

B = [Xd; 
    Md; 
    0];

% Output: u, q, theta, ax
C = [ 0, 1, 0; 
    Xu, Xq, 0]; 

D = [0;
    Xd];


end