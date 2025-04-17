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

function [u_mpc, exitflag] = interiorPoint(H, f, A, b, A_i, b_i, n_u, N)

exitflag = 0;

% initial guesses for Lagrange multipliers (y, z) and slacks (s)
alpha = 1.e-3; % step size
y = ones(size(A, 1), 1);
z = ones(size(A_i, 1), 1);
Z = diag(z);
m = length(b_i);
mu = 100;
sigma = 0.9;
tau = 0.995;

% Phase 1 method for generating feasible initial guess
nx = size(H,1); ns = length(b_i);
f_init = [ones(ns, 1); zeros(nx, 1)];
A_init = [-eye(ns) zeros(ns,nx);
    eye(ns) A_i];
b_init = [-5.*ones(ns,1); b_i];
options = optimoptions('linprog', 'Display', 'off');
init_guess = linprog(f_init,A_init,b_init,[],[],[],[],options);
s = init_guess(1:ns);
S = diag(s);
u = init_guess(ns+1:end);

% anonymous functions
tol = 1.e-12;
c_e = @(x) A*x - b;
c_i = @(x) A_i*x - b_i;
E = @(u, s, z) max([norm(H*u-(A_i'*z)+f, inf) ...
    norm(diag(s)*z - mu.*ones(length(s),1), inf) ...
    norm(c_i(u) + s,inf)]);

% KKT matrices based on existence of equality/inequality constraints
H_kkt = [];
grad_kkt = [];
current_optimum = [];
x0 = [init_guess; z];

while mu > 0.1
    % H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i';
    %     zeros(size(Z,1), size(H,2)) Z S;
    %     A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1))];
    % grad_kkt = [H*u-(A_i'*z)+f; S*z - mu.*ones(size(S,1),1); c_i(u) + s];

    H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i';
        A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1));
        zeros(size(Z,1), size(H,2)) Z S];
    grad_kkt = [H*u-(A_i'*z)+f; c_i(u) + s; S*z - mu.*ones(size(S,1),1)];

    % H_kkt = [A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1));
    %     H zeros(size(H,1), size(Z, 2)) A_i';
    %         zeros(size(Z,1), size(H,2)) Z S];
    % grad_kkt = [c_i(u) + s; H*u-(A_i'*z)+f; S*z - mu.*ones(size(S,1),1)];

    c_i = @(x) A_i*x - b_i;
    E = @(u, s, z) max([norm(H*u-(A_i'*z)+f, inf) ...
        norm(diag(s)*z - mu.*ones(length(s),1), inf) ...
        norm(c_i(u) + s,inf)]);
    while (E(u, s, z) > mu)

        c_i = @(x) A_i*x - b_i;
        E = @(u, s, z) max([norm(H*u-(A_i'*z)+f, inf) ...
            norm(diag(s)*z - mu.*ones(length(s),1), inf) ...
            norm(c_i(u) + s,inf)]);

        % [row, col, val] = find(H_kkt);
        % coo = sortrows([row col val], 1);
        % row = coo(:,1); col = coo(:,2); val = coo(:,3);

        % Improve this
        % p = -H_kkt\grad_kkt;

        %% Gaussian elimination + back-substitution

        % [U,c,~] = GaussElimPivoting(H_kkt,-grad_kkt); % unoptimized
        % [p,~] = backSubs(U,c);

        [U,c,~] = GaussElim_IP(H_kkt,-grad_kkt,N,m); % optimized
        [p,~] = backSubs_IP(U,c,N,m);

        % inv_H_kkt = schurInverse(H_kkt,N,m);
        % inv_H_kkt = schurInverseNewton(H_kkt,N,m);
        % % inv_H_kkt = mtx_inv_newton(H_kkt,1.e-9,100,H_kkt'/(norm(H_kkt,1)*norm(H_kkt,inf)));
        % p = -inv_H_kkt*grad_kkt;
        % p = -H_kkt\grad_kkt;

        %% Iterative methods
        %
        % D = diag(diag(H_kkt));
        % L = -1.*tril(H_kkt,-1);
        % U = -1.*triu(H_kkt, 1);
        % H_SOR = inv(D-w.*L)*((1-w).*D+w.*U);
        % P_cond = diag(1 ./ diag(H_kkt));
        % H_cond = P_cond*H_kkt;
        % D = diag(diag(H_cond));
        % L = -1.*tril(H_cond,-1);
        % U = -1.*triu(H_cond, 1);
        % H_SOR_cond = inv(D-w.*L)*((1-w).*D+w.*U);

        w = 0.54;
        % [p,~,~] = SOR(H_kkt, -grad_kkt, zeros(length(grad_kkt),1), w, 1.e-5, 100);
        % [p,~,~] = SOR_IP(H_kkt, -grad_kkt, zeros(length(grad_kkt),1), w, 1.e-2, 100, N, m);
        % [p,~,~] = SOR_COO(row, col, val, -grad_kkt, zeros(length(grad_kkt),1), w, 1.e-2, 100);

        lu = length(u); ls = length(s); lz = length(z);
        p_u = p(1:lu);
        p_s = p(lu+1:lu+ls);
        p_z = p(lu+ls+1:lu+ls+lz);
        % p_s = p(1:ls);
        % p_u = p(ls+1:ls+lu);
        % p_z = p(ls+lu+1:lu+ls+lz);
        alpha_s = LineSearch(s, p_s, tau);
        alpha_z = LineSearch(z, p_z, tau);
        u = u + alpha_s*p_u;
        s = s + alpha_s*p_s;
        z = z + alpha_z*p_z;
        S = diag(s); Z = diag(z);

        % % update KKT Hessian inverse via Sherman-Woodbury formula
        % H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i';
        %     zeros(size(Z,1), size(H,2)) Z S;
        %     A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1))];
        % grad_kkt = [H*u-(A_i'*z)+f; S*z - mu.*ones(size(S,1),1); c_i(u) + s];

        H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i';
            A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1));
            zeros(size(Z,1), size(H,2)) Z S];

        grad_kkt = [H*u-(A_i'*z)+f; c_i(u) + s; S*z - mu.*ones(size(S,1),1)];

        %     H_kkt = [A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1));
        %     H zeros(size(H,1), size(Z, 2)) A_i';
        %         zeros(size(Z,1), size(H,2)) Z S];
        % grad_kkt = [c_i(u) + s; H*u-(A_i'*z)+f; S*z - mu.*ones(size(S,1),1)];
    end
    mu = sigma*mu; % tighten mu
end
u_mpc = u(1:n_u);
end

