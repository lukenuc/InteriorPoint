
% Test script for the SOR function optimized for interior point. The result
% here indicates that w = 1 (Gauss-Seidel) converges the quickest. 

clear; clc; close all; 
include
load('test.mat')

w_arr = 0:0.01:2; 
times = zeros(1, length(w_arr)); 
iters = zeros(1, length(w_arr));
num_trials = 50; 

for i = 1:length(w_arr)
    i
    time = []; 
    for j = 1:num_trials
        tic
        [p,k,~] = SOR_IP(H_kkt, -grad_kkt, zeros(length(grad_kkt),1), w_arr(i), 1.e-9, 100, N, m);
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