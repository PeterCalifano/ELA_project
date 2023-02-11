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

% addpath('datasets','common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');

DefaultFontSize = 16;
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaulttextinterpreter', 'latex');
set(0, 'defaultAxesFontSize', DefaultFontSize)

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

% Delays
% Delay = 4 sampling intervals = 4 ms
delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

%% Load controller parameters

parameters_controller                    

%% M injection example (sweeep: first column time vector, second column time history of pitching moment) 
load ExcitationM
% SetPoint for Translation control
SetPoint = [0, 0];

%% Values selected
% Time grid
t = ExcitationM(:, 1);
simulation_time = t(end) - t(1);
decimation = 1; % [s]

%% Launch SIMULATOR
save('input_workspace.mat');
model_name = 'Simulator_Single_Axis';

% Simulate or load sample output
% if exist('simout.mat', 'file')
output = sim(model_name, 'TimeOut', simulation_time);
save('simout.mat', "output");
% else
%     load('simout.mat');
% end

%% Delete temporary files

if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end

%% Signals Pre-Processing
Excit_signal = output.Excit_signal;
RemoveZeroInputMask = Excit_signal ~= 0;

% Extract useful input/output samples 
Mtot = output.Mtot(RemoveZeroInputMask);
time_grid = output.time_grid(RemoveZeroInputMask);
ax = output.ax(RemoveZeroInputMask);
q = output.q(RemoveZeroInputMask);

% dt = 1/250; % 250 Hz, defined in parameters_controller
time_grid = time_grid(5:end);
% Consider delay of the output (4 samples)
Mtot = Mtot(1:end-4);
ax = ax(5:end);
q = q(5:end);

% TO MODIFY: DELAY APPLIES TO TIME GRID, not SAMPLE INDEX. Must be
% equivalent to: 

% s_time_out = s_time + 1/fs*4; %sampling time for the output + delay
% % Normalized pitch moment resampling at 250 Hz
% M = M - mean(M); %removed bias
% meas.M = interp1(tM,M,s_time_in);

% Remove time samples where Excitation signal 

% figure;
% plot(time_grid, Mtot, '-');
% xlabel('Time [s]');
% ylabel('Excitation signal')
% grid minor;
% 
% figure;
% plot(time_grid, ax, '-');
% xlabel('Time [s]');
% ylabel('Longitudinal acceleration [m/s^2]')
% grid minor;
% 
% figure;
% plot(time_grid, rad2deg(q), '-');
% xlabel('Time [s]')
% ylabel('Pith rate q [deg/s]');
% grid minor;
% 
% close all

% TO DO: EVALUATE WHETHER LOWPASS FILTERING IS USEFUL
% Call script to estimate Frequency Response Function from output time signals
% output_delay = 0.016; % [s]

% Call script for Model Identification
IdentifyModel;



%% END OF CODE





