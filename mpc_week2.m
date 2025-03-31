% Week 1 (Jan 28 - Feb 3): Constrained MPC for Tracking
% 
% Objective: Extend the controller to handle reference tracking with input and state constraints.
% Tasks:
% Modify your system to track a sinusoidal reference signal.
% Add input/output constraints (e.g., actuator limits).
% Deliverable: Demonstrate tracking performance under constrained conditions.

clear all
close all
clc
x = zeros(2,1); % state vector 
u = 0; % control input 
Ts = 0.1; % sampling period [s] 
Tsim = 20; % simulation period [s]
t = 0:Ts:Tsim;

b = 0.05; % drag coefficient
A = [0 1; 0 -b]; % Continuous time state-space matrices for double integrator
B = [0; 1]; 
C = eye(2);
D = [0; 0]; 

Ad = eye(2) + Ts.*A; % Discrete time dynamics derived via forward Euler integration
Bd = Ts.*B; 

%% MPC controller using MPC toolbox 

plant = ss(A, B, C, D); 
% mpcDesigner % optional visual tool
p = 10; % prediction horizon
m = 5; % control horizon
mpcController = mpc(plant, Ts, p, m); % create the MPC controller instance

% Output and manipulated variable constraints
mpcController.OV(1).Min = -Inf;
mpcController.OV(1).Max = 100; % constraint on position [m], need car to park here
mpcController.OV(2).Min = -Inf;  
mpcController.OV(2).Max = 22; % constraint on velocity [m/s], can't go too fast
mpcController.MV(1).Min = -10; 
mpcController.MV(1).Max = 10; % constraint on applied force [m/s^2] for comfortable driving experience

% Adjusting weights for each variable
mpcController.Weights.OV = [100, 10]; % [w_pos, w_vel]
mpcController.Weights.MV = 10; % [w_u]
mpcController.Weights.MVRate = 1; % [w_delta_u]

% Controller simulation
state = mpcstate(mpcController); 
ref = [100; 0]; 
mpc_states = zeros(2, length(t)); 
u_mpc = zeros(1, length(t)); 
mpc_states(:, 1) = [0; 0]; 
x = [0; 0];

t_bez = linspace(0, 1, Tsim/Ts);
bezier_ref = kron((1-t_bez).^3, 0) + kron(3*(1-t_bez).^2.*t_bez, -5) + kron(3*(1-t_bez).*t_bez.^2, 105) + kron(t_bez.^3, 100); 
% plot(t, bezier_ref)


for i = 1:length(t)-1
    y = C*x; % measure the outputs
    u = mpcmove(mpcController, state, y, [bezier_ref(i); 0]); % get the MPC control input
    x = Ad*x + Bd*u; % this should be discrete-time dynamics
    mpc_states(:, i+1) = x; 
    u_mpc(i) = u; 
end

% plot results 
figure('Position', [100 100 1200 400]); 
subplot(1, 3, 1)
stairs(t(1:end-1), u_mpc(1:end-1), 'LineWidth', 2); 
title('Vehicle Control: MPC Input'); 
xlabel('Time [s]')
ylabel('Applied Acceleration [m/s/s]')
grid on 
grid minor
subplot(1, 3, 2)
plot(t(1:end-1), mpc_states(1,1:end-1), 'LineWidth', 2, 'Color', 'g'); 
hold on
plot(t(1:end-1), bezier_ref, 'b--')
hold off
title('Vehicle Control: Position'); 
xlabel('Time [s]')
ylabel('Position [m]')
grid on 
grid minor
subplot(1, 3, 3)
plot(t(1:end-1), mpc_states(2,1:end-1), 'LineWidth', 2, 'Color', 'r'); 
title('Vehicle Control: Velocity'); 
xlabel('Time [s]')
ylabel('Velocity [m/s]')
grid on 
grid minor

%% Car visualization

figure; 
hold on;
xlim([-10, 110]); 
ylim([-60, 60]); 
car = plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Car representation
car_width = 4.9; car_height = 1.94; 
for i = 1:length(t(1:end-1))
    cla;
    rectangle('Position', [mpc_states(1,i), car_height / 2, car_width, car_height], ...
              'Curvature', 0.2, ...
              'FaceColor', 'r');
    grid on
    title('MPC Control Car Simulation')
    xlabel('Horizontal Distance [m]')
    ylabel('Vertical Distance [m]')
    pause(Ts); % Adjust the pause for animation speed
end

hold off;