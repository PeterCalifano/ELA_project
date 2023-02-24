function [Mtot, ax, q, time_grid] = OutputPreProcess(output, N_delay)

% Extract useful input/output samples
Mtot = output.Mtot;
time_grid = output.time_grid;
ax = output.ax;
q = output.q;

% dt = 1/250; % 250 Hz, defined in parameters_controller
time_grid = time_grid((1+N_delay):end);
% Consider delay of the output (4 samples)
Mtot = Mtot(1:(end-N_delay));
ax = ax((1+N_delay):end);
q = q((1+N_delay):end);

end