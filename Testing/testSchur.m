
% Test script for the Schur complement for computing inverses

clear all; clc; close all;
include
load('test.mat')

w_arr = 1;
times = zeros(1, length(w_arr));
iters = zeros(1, length(w_arr));
num_trials = 50;

for i = 1:length(w_arr)
    i
    time = [];
    for j = 1:num_trials
        tic
        inv_H_kkt = schurWoodbury(H_kkt,N,m);
        % inv_H_kkt = schurInverse(H_kkt,N,m);
        p = -inv_H_kkt*grad_kkt;
        time = [time toc];
    end
    time(time == max(time) | time == min(time)) = [];
    times = mean(time)
end
% plot(w_arr, times)
% xlabel('\omega'); ylabel('Avg. Time [s]');
% figure
% plot(w_arr, iters)
% xlabel('\omega'); ylabel('Iterations')