
% Test script for the SOR function optimized for interior point. The result
% here indicates that w = 0.4-0.5 (an underrelaxed system) converges the
% quickest. 

clear; clc; close all;
include
load('test.mat')

w_arr = 0:0.01:2;
times = zeros(1, length(w_arr));
iters = zeros(1, length(w_arr));
num_trials = 100;
Sigma = mu.*(inv(S))^2;
H_sym = [H + 2*A_i'*Sigma*A_i -A_i'; -A_i inv(Sigma)];
grad_kkt = [H*u-(A_i'*z)+f; c_i(u) - mu*inv(Z)*ones(m,1)];
r1 = grad_kkt(1:N); r2 = grad_kkt(N+1:end); 
grad_kkt = [r1 - 2*A_i'*Sigma*r2; r2]; 
n = length(grad_kkt); 
[row, col, v] = find(H_sym);
coo = sortrows([row col v], 1);
row = coo(:,1); col = coo(:,2); v = coo(:,3);

for i = 1:length(w_arr)
    i
    time = [];
    for j = 1:num_trials
        tic
        % if j == 1
            [p,k,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), w_arr(i), 1.e-3, 500);
        % else
        %     [p,k,~] = SOR_COO(row, col, v, -grad_kkt, p + rand(n,1), w_arr(i), 1.e-3, 1000);
        % end
        % [p,k,~] = SOR_IP(H_kkt, -grad_kkt, zeros(length(grad_kkt),1), w_arr(i), 1.e-3, 1000, N, m);
        % [p,k,~] = SOR(H_kkt,-grad_kkt,zeros(length(grad_kkt),1),1,1.e-5,1000);
        time = [time toc];
    end
    time(time == max(time) | time == min(time)) = [];
    times(i) = mean(time);
    iters(i) = k;
end
plot(w_arr, times)
xlabel('\omega'); ylabel('Avg. Time [s]');
figure
plot(w_arr, iters)
xlabel('\omega'); ylabel('Iterations')