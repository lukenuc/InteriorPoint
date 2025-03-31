% Week 1 (Jan 20-27): Linear MPC Basics
% 
% Objective: Implement and simulate a basic MPC controller for a simple linear system (e.g., double integrator or first-order plant).
% Tasks:
% Write state-space models in MATLAB.
% Design and simulate an MPC controller using MATLAB’s Model Predictive Control Toolbox.
% Explore constraints on states and inputs.
% Deliverable: MATLAB script or Simulink model with plots showing system response with and without MPC.

clear all
close all
clc
x = zeros(2,1); % state vector 
u = 0; % control input 
Ts = 0.1; % sampling period [s] 
Tsim = 20; % simulation period [s]

b = 0.05; % drag coefficient
A = [0 1; 0 -b]; % Continuous time state-space matrices for double integrator
B = [0; 1]; 
C = eye(2);
D = [0; 0]; 

Ad = eye(2) + Ts.*A; % Discrete time dynamics derived via forward Euler integration
Bd = Ts.*B; 

%% Demonstrating the zero-input and zero-state response of our system

t = 0:Ts:Tsim; 
zsr_states = zeros(2, length(t)); 
zir_states = zeros(2, length(t)); 

% zero-input response
u = 0;
x = [100; 5]; % Start at x = 100 [m] traveling at v = 5 [m/s] 
zir_states(:, 1) = x;
for i = 1:length(t)-1
    x = Ad*x + Bd*u; 
    zir_states(:, i+1) = x; 
end

% zero-state response
u = -0.5;
x = [0; 0]; % Deccelerate at -5 [m/s^2] 
zsr_states(:, 1) = x;
for i = 1:length(t)-1
    x = Ad*x + Bd*u; 
    zsr_states(:, i+1) = x; 
end

% plot results 
figure('Position', [100 100 1000 500]); 
subplot(1, 2, 1)
plot(t, zir_states(1,:), 'LineWidth', 2, 'Color', 'g'); 
title('Double Integrator: Zero-Input Response'); 
xlabel('Time [s]')
ylabel('Position [m]')
grid on 
subplot(1, 2, 2)
plot(t, zsr_states(1,:), 'LineWidth', 2, 'Color', 'g'); 
title('Double Integrator: Zero-State Response'); 
xlabel('Time [s]')
ylabel('Position [m]')
grid on 

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
mpcController.Weights.OV = [2, 2]; % [w_pos, w_vel]
mpcController.Weights.MV = 10; % [w_u]
mpcController.Weights.MVRate = 0.1; % [w_delta_u]

% Controller simulation
state = mpcstate(mpcController); 
ref = [100; 0]; 
mpc_states = zeros(2, length(t)); 
u_mpc = zeros(1, length(t)); 
mpc_states(:, 1) = [0; 0]; 
for i = 1:length(t)-1
    y = C*x; % measure the outputs
    u = mpcmove(mpcController, state, y, ref); % get the MPC control input
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
xlim([-120, 120]); 
ylim([-120, 120]); 
car = plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Car representation
car_width = 20; car_height = 8; 
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