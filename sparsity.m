clear all; close all; clc;

Ts = 0.1; % sampling period [s]
b = 0.05; % drag coefficient
A = [0 1; 0 -b];
B = [0; 1];
C = eye(2);
D = [0; 0];
Ad = eye(2) + Ts.*A;
Bd = Ts.*B;
nx = 2; % number of states
nu = 1; % number of control inputs
x_ref = [10; 0];

p_values = 1:20; % range of p values to test
sparsity_percent = zeros(size(p_values));
selected_p = [2, 5, 10, 20]; % Chosen p values for sparsity map

figure;
tiledlayout(2,2);

for idx = 1:length(p_values)
    p = p_values(idx);
    x0_t = [0; 0];
    xmin = repmat([-100; -30], p, 1);
    xmax = repmat([100; 22], p, 1);
    umin = repmat([-10], p, 1);
    umax = repmat([10], p, 1);
    
    qy_weight = [50; 1];
    qu_weight = 1;
    Qy = diag(qy_weight);
    Qu = diag(qu_weight);
    My = diag(repmat(qy_weight, p, 1));
    Mu = diag(repmat(qu_weight, p, 1));
    Mc = kron(eye(p), C);
    
    M_ab = zeros(nx*p, nu*p);
    M_ak = zeros(nx*p, nx);
    
    for i = 1:p
        M_ak((i-1)*nx+1:i*nx, :) = C*Ad^i;
        for j = 1:p
            if j > i
                M_ab((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = zeros(nx, nu);
            else
                M_ab((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = C*(Ad^(i-j))*Bd;
            end
        end
    end
    
    A_i = [M_ab; -M_ab; eye(nu*p); -eye(nu*p)];
    b_i = [xmax-M_ak*(x0_t-x_ref); -xmin+M_ak*(x0_t-x_ref); umax; -umin];
    H = (1/2) * ((M_ab')*My*M_ab + Mu);
    
    z = ones(size(A_i, 1), 1);
    Z = diag(z);
    s = ones(size(b_i, 1), 1);
    S = diag(s);
    H_kkt = [H zeros(size(H,1), size(Z, 2)) A_i';  
             zeros(size(Z,1), size(H,2)) Z S; 
             A_i eye(size(A_i, 1), size(Z, 2)) zeros(size(A_i,1), size(A_i,1))];
    
    sparsity_percent(idx) = 100 * (1 - nnz(H_kkt) / numel(H_kkt));
    
    if ismember(p, selected_p)
        nexttile;
        spy(H_kkt);
        title(sprintf('Sparsity for p = %d', p));
    end
end

%%
figure;
plot(p_values, sparsity_percent, 'o-', 'LineWidth', 2, 'Color', 'b');
xlabel('Prediction Horizon (p)');
ylabel('Sparsity (%)');
title('Sparsity of H_{kkt} vs. Prediction Horizon');
grid on;
