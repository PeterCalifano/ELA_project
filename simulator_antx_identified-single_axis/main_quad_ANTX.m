%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANT-X SIMULATOR - MAIN                                                  %
% Authors:  Mattia Giurato (mattia.giurato@polimi.it)                     %
%           Paolo Gattazzo (paolo.gattazzo@polimi.it)                     %
% Date: 13/12/2017                                                        %
% Adapted to ANT-X 2DoF by:  Salvatore Meraglia (salvatore.meraglia@polimi.it)%
% Date: 22/12/2022                                                        %
%
% Further modified to include structure three-state identified longitudinal model
% 06/01/23 ML
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clearvars;
close all;
clear;
clc;

addpath('common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');

DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize);

rng default;


%% Model parameters
% Initial model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)
 
% Derivative of udot wrt u
Xu = -0.1068;
% Derivative of udot wrt q
Xq = 0.1192;
% Derivative of qdot wrt u
Mu = -5.9755;
% Derivative of qdot wrt q
Mq = -2.6478;
% Derivative of udot wrt delta
Xd = -10.1647;
% Derivative of qdot wrt delta
Md = 450.71;

th_true = [Xu, Xq, Mu, Mq, Xd, Md];

%% State space model
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

s = tf('s');

% Noise
% noise.Enabler = 0;
noise.Enabler = 1;

noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]

noise.vel_stand_dev = noise.Enabler * 0.01;                                %[m/s]

noise.attitude_stand_dev = noise.Enabler * deg2rad(0.0076);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(0.01);                   %[rad/s]

% Delays (samples)
delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

%% Load controller parameters
parameters_controller                    

%% M injection example (sweeep: first column time vector, second column time history of pitching moment) 
load ExcitationM
% SetPoint for Translation control
SetPoint = [0, 0];

%% Launch SIMULATOR
model_name = 'Simulator_Single_Axis';
save('input_workspace.mat');

% Time grid
t = ExcitationM(:, 1);
simulation_time = t(end) - t(1);
decimation = 1; % [s]

% Simulate system and save output
output = sim(model_name, 'TimeOut', simulation_time);
save('simout.mat', "output");


%% Delete temporary files

if exist('/slprj', 'dir')
    rmdir('slprj', 's')                                                    
end

%% Signals Pre-Processing
N_delay = 1;

time_grid = output.time_grid;
CutInputMask = time_grid >= 23 & time_grid <= 100;

% Extract useful input/output samples 
Excit_signal = output.Excit_signal(CutInputMask);
Mtot = output.Mtot(CutInputMask);
ax = output.ax(CutInputMask);
q = output.q(CutInputMask);
time_grid = time_grid(CutInputMask);

% dt = 1/250; % 250 Hz, defined in parameters_controller
% Consider delay of the output (1 samples): shift of I/O signals
time_grid = time_grid((1+N_delay):end);
Excit_signal = Excit_signal(1:(end-N_delay));
Mtot = Mtot(1:(end-N_delay));
ax = ax((1+N_delay):end);
q = q((1+N_delay):end);

%% Input and Output signals plot

figure;
plot(time_grid, Mtot, '-', 'LineWidth', 1.02);
hold on;
plot(time_grid, Excit_signal, '-', 'LineWidth', 1.05)
xlabel('Time [s]');
ylabel('Signal Amplitude [-]')
% Default options
grid minor
axis auto
ax_gca = gca;
ax_gca.XAxisLocation = 'bottom';
ax_gca.YAxisLocation = 'left';
ax_gca.XMinorTick = 'on';
ax_gca.YMinorTick = 'on';
ax_gca.LineWidth = 1.04;
hold off;
legend('Excitation signal', 'Input torque')

figure;
plot(time_grid, ax, '-', 'LineWidth', 1.02);
xlabel('Time [s]');
ylabel('Longitudinal acceleration [$m/s^2$]')
% Default Options
grid minor
axis auto
ax_gca = gca;
ax_gca.XAxisLocation = 'bottom';
ax_gca.YAxisLocation = 'left';
ax_gca.XMinorTick = 'on';
ax_gca.YMinorTick = 'on';
ax_gca.LineWidth = 1.04;
hold off;

figure;
plot(time_grid, rad2deg(q), '-', 'LineWidth', 1.02);
xlabel('Time [s]')
ylabel('Pith rate q [deg/s]');
% Default Options
grid minor
axis auto
ax_gca = gca;
ax_gca.XAxisLocation = 'bottom';
ax_gca.YAxisLocation = 'left';
ax_gca.XMinorTick = 'on';
ax_gca.YMinorTick = 'on';
ax_gca.LineWidth = 1.04;
hold off;

% Call script for Model Identification
run('IdentifyModel');
clear;

return
% Call script for Optimization of the input signal
run('OptimizeIdentification.m');

%% END OF CODE





