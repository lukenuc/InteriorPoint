% =========================================================================
% FILE NAME:    profiling.m
% AUTHOR:       Luke Nuculaj
%
% DESCRIPTION: This script collects a timing profile of several versions of
% the "interiorPoint" function, and uses MATLAB's "quadprog" function as a
% timing benchmark. As a proof-of-concept, a position reference tracking
% problem with simplified vehicle dynamics is utilized.
% =========================================================================

clear all
close all
clc

p_arr = 1:1:30;
avg_times = zeros(4, length(p_arr));

for count = 1:length(p_arr)
    count
    x = zeros(2,1); % state vector
    u = 0; % control input
    Ts = 0.1; % sampling period [s]
    Tsim = 0.5; % simulation period [s]
    t = 0:Ts:Tsim;

    b = 0.05; % drag coefficient
    A = [0 1; 0 -b]; % Continuous time state-space matrices for double integrator
    B = [0; 1];
    C = eye(2);
    D = [0; 0];

    Ad = eye(2) + Ts.*A; % Discrete time dynamics derived via forward Euler integration
    Bd = Ts.*B;

    %% MPC setup

    nx = 2; % number of states
    nu = 1; % number of control inputs
    p = p_arr(count); % size of prediction horizon
    x0_t = [0; 0]; % initial state

    % x0_t = repmat(x0, p, 1);
    xmin = repmat([-100; -30], p, 1);
    xmax = repmat([100; 22], p, 1);
    umin = repmat([-10], p, 1);
    umax = repmat([10], p, 1);
    M_ab = zeros(nx*p, nu*p);
    M_ak = zeros(nx*p, nx);
    qy_weight = [50; 1];
    qu_weight = 1;
    Qy = diag(qy_weight);
    Qu = diag(qu_weight);
    My = diag(repmat(qy_weight, p, 1));
    Mu = diag(repmat(qu_weight, p, 1));
    Mc = kron(eye(p), C);

    % assembling the dense matrices
    for i = 1:p
        M_ak((i-1)*nx+1:(i-1)*nx+nx, 1:nx) = C*Ad^i;
        for j = 1:p
            if (j > i)
                M_ab((i-1)*nx+1:(i-1)*nx+nx, (j-1)*nu+1:(j-1)*nu+nu) = zeros(nx,nu);
            else
                M_ab((i-1)*nx+1:(i-1)*nx+nx, (j-1)*nu+1:(j-1)*nu+nu) = C*(Ad^(i-j))*Bd;
            end
        end
    end

    % computing matrices for constraints and objective function
    x_ref = [10; 0]; % desired location

    A_i = [M_ab; -M_ab; eye(nu*p); -eye(nu*p)];
    b_i = [xmax-M_ak*(x0_t-x_ref); -xmin+M_ak*(x0_t-x_ref); umax; -umin];
    H = (1/2).*((M_ab')*My*M_ab + Mu);
    % f = (M_ab')*My*M_ak*(x0_t-x_ref);

    %% MPC controller using MPC toolbox

    % Reference Bezier curve
    t_bez = linspace(0, 1, Tsim/Ts);
    t_bez = [t_bez repmat(t_bez(end), 1, p)]; % extend ending for prediction horizon
    bezier_ref = kron((1-t_bez).^3, 0) + kron(3*(1-t_bez).^2.*t_bez, -5) + kron(3*(1-t_bez).*t_bez.^2, 105) + kron(t_bez.^3, 100);

    for i = 1:length(t)-1

        pos_ref = bezier_ref(i:i+p-1)';
        y_ref = zeros(nx*p, 1);
        y_ref(1:2:end) = pos_ref;
        f = (M_ab')*My*(M_ak*x0_t - y_ref);
        b_i = [xmax-M_ak*x0_t+y_ref; -xmin+M_ak*x0_t-y_ref; umax; -umin];

        % Naive; inverts KKT Hessian each time
        tic
        [u,~] = interiorPointInv(H,f,[],[],A_i,b_i,nu);
        avg_times(1, count) = avg_times(1, count) + toc/(length(t)-1);

        % % Gaussian Elim; solves linear system with GE + BS each step
        % tic
        % [u,~] = interiorPointGE(H,f,[],[],A_i,b_i,nu);
        % avg_times(2, count) = avg_times(2, count) + toc/(length(t)-1);

        % Sherman-Woodbury; inverts KKT Hessian once and does rank one
        % updates
        tic
        [u,~] = interiorPoint(H,f,[],[],A_i,b_i,nu);
        avg_times(3, count) = avg_times(3, count) + toc/(length(t)-1);

        % quadprog; MATLAB benchmark
        opts = optimoptions(@quadprog,...
            'Algorithm','interior-point-convex');
        tic
        [u,~] = quadprog(H,f,A_i,b_i,[],[],[],[],x0_t,opts);
        avg_times(4, count) = avg_times(4, count) + toc/(length(t)-1);

        % progressBar(floor(100*i/(length(t)-1)))
        % pos_ref = bezier_ref(i:i+p-1)';
        % y_ref = zeros(nx*p, 1);
        % y_ref(1:2:end) = pos_ref;
        % f = (M_ab')*My*(M_ak*x0_t - y_ref);
        % b_i = [xmax-M_ak*x0_t+y_ref; -xmin+M_ak*x0_t-y_ref; umax; -umin];
        % [u,~] = interiorPoint(H,f,[],[],A_i,b_i,nu);
        % x0_t = Ad*x0_t + Bd*u; % update state with optimal control
        %
        % mpc_states(:, i+1) = x0_t;
        % u_mpc(i) = u;
    end
end