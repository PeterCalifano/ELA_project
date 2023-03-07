function sim_object = SetModel(theta, ExcitationM, noise_flag)

% Function to set the Simulation object to execute simulations
model_name = 'Simulator_Single_Axis';
% Time grid
t = ExcitationM(:, 1);
simulation_time = t(end) - t(1);

if exist('noise_flag', 'var')
    noise.Enabler = noise_flag;
    noise.pos_stand_dev = noise.Enabler * 0.0011;                            	% [m]
    noise.vel_stand_dev = noise.Enabler * 0.01;                                 % [m/s]
    noise.attitude_stand_dev = noise.Enabler * deg2rad(0.0076);                 % [rad]
    noise.ang_rate_stand_dev = noise.Enabler * deg2rad(0.01);                   % [rad/s]
else
    % Do nothing
    noise.Enabler = 1;
    noise.pos_stand_dev = noise.Enabler * 0.0011;                            	% [m]
    noise.vel_stand_dev = noise.Enabler * 0.01;                                 % [m/s]
    noise.attitude_stand_dev = noise.Enabler * deg2rad(0.0076);                 % [rad]
    noise.ang_rate_stand_dev = noise.Enabler * deg2rad(0.01);                   % [rad/s]
end

% Determine State space model
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
sim_object = setVariable(sim_object, 'noise', noise);

end