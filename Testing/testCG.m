
% Test script for the Schur complement for computing inverses

clear all; clc; close all;
include
load('test.mat')

times = 0;
num_trials = 2000;

Sigma = mu.*(inv(S))^2;
H_sym = [H + 2*A_i'*Sigma*A_i -A_i'; -A_i inv(Sigma)];
grad_kkt = [H*u-(A_i'*z)+f; c_i(u) - mu*inv(Z)*ones(m,1)];
r1 = grad_kkt(1:N); r2 = grad_kkt(N+1:end); 
grad_kkt = [r1 - 2*A_i'*Sigma*r2; r2]; 
% M = 1 ./ sqrt(diag(H_sym)); % preconditioner
M = ones(size(H_sym,1),1); 
[row, col, v] = find(H_sym);
coo = sortrows([row col v], 1);
row = coo(:,1); col = coo(:,2); v = coo(:,3);
n = length(grad_kkt); 

time = [];
for j = 1:num_trials
    tic
    % [p, k] = conjugateGradient(H_sym, -grad_kkt, M, zeros(length(grad_kkt), 1), 1e-5, 100);
    [p,k,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), 0.45, 1.e-9, 10);
    time = [time toc];
end
time(time == max(time) | time == min(time)) = [];
times = mean(time)
p_temp = -H_sym\grad_kkt; 
norm(p-p_temp,inf)