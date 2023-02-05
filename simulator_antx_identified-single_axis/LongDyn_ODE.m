function [A, B, C, D] = LongDyn_ODE(Xu, Xq, Mu, Mq, Xd, Md, Ts)

% Initial model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)

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