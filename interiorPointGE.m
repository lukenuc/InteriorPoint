% =========================================================================
% FUNCTION NAME: interiorPoint.m
% AUTHOR:        Luke Nuculaj
% DESCRIPTION:   Takes the matrices for a finite-horizon OCP problem of
% interest as input and performs several interior point iterations,
% tightening the relaxation parameter \mu with each step. Returns the
% approximate optimal control input arrived at. 
%
% INPUTS:
%   - H : objective function Hessian matrix
%   - f : objective function gradient
%   - A: equality constraint gradient
%   - b: equality constraint bound
%   - A_i: inequality constraint gradient
%   - b_i: inequality constraint bound
%
% OUTPUTS: 
%   - u_mpc: optimal control move
%
% EXIT FLAGS: 
%   - exitflag: = 0 means the function was successful
% 
% =========================================================================

function [u_mpc, exitflag] = interiorPoint(H, f, A, b, A_i, b_i, n_u)
    
    exitflag = 0; 
    
    % convergence criteria (check infinity norm on complementary slackness)
    tol = 1.e-3; 
    is_converged = @(Sz, tol) (norm(Sz, inf) <= tol); 
    
    % anonymous functions for conveniently evaluating constraints
    c_e = @(x) A*x - b; 
    c_i = @(x) A_i*x - b_i; 

    % initial guesses for Lagrange multipliers (y, z) and slacks (s)
    alpha = 1.e-3; % step size
    u = ones(size(H,1), 1); 
    y = ones(size(A, 1), 1); 
    z = ones(size(A_i, 1), 1); 
    Z = diag(z); 
    s = ones(size(b_i, 1), 1); 
    S = diag(s); 
    mu = 10; 
    % KKT matrices based on existence of equality/inequality constraints
    H_kkt = []; 
    grad_kkt = []; 
    current_optimum = []; 
    if (isempty(A))
        H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i'; 
                zeros(size(Z,1), size(H,2)) Z S; 
                A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1))];
        % grad_kkt = [zeros(size(H,1),1); S*z - mu*ones(size(S,1),1); c_i(u) + s];
        grad_kkt = [H*u-(A_i'*z)+f; S*z - mu.*ones(size(S,1),1); c_i(u) + s];
        current_optimum = [u; s; z]; 
    elseif (isempty(A_i))

    else

    end
    
    % inv_H_kkt = inv(H_kkt); 
    % U = [zeros(size(H,1),size(Z,1)); eye(size(Z,1)); zeros(size(A_i, 1), size(Z,1))];
    while (~is_converged(S*z, tol))
        [H_kkt, grad_kkt] = GaussElimPivoting(H_kkt, grad_kkt); 
        [p, ~] = backSubs(H_kkt, grad_kkt); 
        % p = inv_H_kkt*grad_kkt; 
        current_optimum = current_optimum - alpha.*p; 
        lu = length(u); ls = length(s); lz = length(z); 
        u = current_optimum(1:lu); 
        s = current_optimum(lu+1:lu+ls); % ensure positivity
        z = current_optimum(lu+ls+1:lu+ls+lz); % ensure positivity
        % delta_S = diag(s) - S; 
        % delta_Z = diag(z) - Z; 
        S = diag(s); Z = diag(z); 
        % update KKT Hessian inverse via Sherman-Woodbury formula
        % V = [zeros(size(Z,1), size(H,2)) delta_Z delta_S];
        % inv_H_kkt = inv_H_kkt*(eye(size(inv_H_kkt,1)) - U*((eye(size(Z,1)) + V*inv_H_kkt*U)\(V*inv_H_kkt)));

        H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i'; 
                zeros(size(Z,1), size(H,2)) Z S; 
                A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1))];
        grad_kkt = [H*u-(A_i'*z)+f; S*z - mu.*ones(size(S,1),1); c_i(u) + s];
        mu = 0.5*mu; % tighten mu 
    end
    u_mpc = u(1:n_u);
end

