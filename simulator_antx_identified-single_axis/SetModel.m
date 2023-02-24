function sim_object = SetModel(theta, ExcitationM)

model_name = 'Simulator_Single_Axis';
% Time grid
t = ExcitationM(:, 1);
simulation_time = t(end) - t(1);

% Simulate system

% [A, B, C, D] = LongDyn_ODE(theta(1), theta(2), theta(3), theta(4), theta(5), theta(6));

Xu = theta(1);
Mu = theta(2);
Xq = theta(3);
Mq = theta(4);
Xd = theta(5);
Md = theta(6);

A = [Xu, Xq, -9.81; 
    Mu, Mq, 0; 
    0, 1, 0];

B = [Xd; 
    Md; 
    0];

% Output: u, q, theta, ax
C = [1, 0, 0; 
    0, 1, 0; 
    0, 0, 1; 
    Xu, Xq, 0]; 

D = [0; 
    0;
    0; 
    Xd];

sim_object = Simulink.SimulationInput(model_name);
sim_object = setVariable(sim_object, 'ExcitationM', ExcitationM);
sim_object = sim_object.setModelParameter('StopTime', num2str(simulation_time));

sim_object = setVariable(sim_object, 'A', A);
sim_object = setVariable(sim_object, 'B', B);
sim_object = setVariable(sim_object, 'C', C);
sim_object = setVariable(sim_object, 'D', D);


end