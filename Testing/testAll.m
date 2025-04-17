
% Test script for the Schur complement for computing inverses

clear all; clc; close all;
include
load('test.mat')

%% Testing execution time on problem of a fixed size (p = 20)

num_methods = 4;
num_trials = 2000;
times = zeros(num_methods, 1);
time_arr = zeros(num_methods, num_trials);
iters = zeros(num_methods, 1);
[row, col, v] = find(H_kkt);
coo = sortrows([row col v], 1);
row = coo(:,1); col = coo(:,2); v = coo(:,3);
p = zeros(length(grad_kkt),1); 
n = length(grad_kkt); 

for j = 1:num_trials

    % Gaussian elimination
    tic
    [U,c,~] = GaussElim_IP(H_kkt,-grad_kkt,N,m);
    [p,~] = backSubs_IP(U,c,N,m);
    time_arr(1, j) = toc;

    % SOR
    tic
    [p,~,~] = SOR_COO(row, col, v, -grad_kkt, p + (1.e-1)*rand(n,1), 0.45, 1.e-9, 100);
    time_arr(2, j) = toc;

    % Schur-Newton
    tic
    inv_H_kkt = schurInverseNewton(H_kkt,N,m);
    p = -inv_H_kkt*grad_kkt;
    time_arr(3, j) = toc;

    % MATLAB built-in inverse
    tic
    p = -H_kkt\grad_kkt;
    time_arr(4, j) = toc;
end

for i = 1:num_methods
    t = time_arr(i,:);
    t(t == max(t) | t == min(t)) = [];
    times(i) = mean(t);
end

% Plot execution times as a bar graph
figure;
methods = {'GE', 'SOR', 'Schur-Newton', 'A\b'};
bar(times);
set(gca, 'XTickLabel', methods, 'XTick', 1:num_methods);
xlabel('Method');
ylabel('Execution Time (seconds)');
title('Comparison of Execution Times for Different Methods');
grid on;

%% Testing execution times of optimized versus unoptimized methods

num_methods = 3;
num_trials = 1000;
times_opt = zeros(num_methods, 1);
times = zeros(num_methods, 1);
time_arr_opt = zeros(num_methods, num_trials);
time_arr = zeros(num_methods, num_trials);
[row, col, v] = find(H_kkt);
coo = sortrows([row col v], 1);
row = coo(:,1); col = coo(:,2); v = coo(:,3);
n = length(grad_kkt); 

for j = 1:num_trials

    % Gaussian elimination
    [U,c,~] = GaussElimPivoting(H_kkt,-grad_kkt); % unoptimized
    [p,~] = backSubs(U,c);
    time_arr(1, j) = toc;
    tic
    [U,c,~] = GaussElim_IP(H_kkt,-grad_kkt,N,m);
    [p,~] = backSubs_IP(U,c,N,m);
    time_arr_opt(1, j) = toc;

    % SOR
    tic
    [p,~,~] = SOR(H_kkt, -grad_kkt, zeros(length(grad_kkt),1), 0.45, 1.e-5, 100);
    time_arr(2, j) = toc;
    tic
    [p,~,~] = SOR_COO(row, col, v, -grad_kkt, zeros(n,1), 0.45, 1.e-5, 100);
    time_arr_opt(2, j) = toc;

    % Schur-Newton
    tic
    inv_H_kkt = schurInverse(H_kkt,N,m);
    p = -inv_H_kkt*grad_kkt;
    time_arr(3, j) = toc;
    tic
    inv_H_kkt = schurInverseNewton(H_kkt,N,m);
    p = -inv_H_kkt*grad_kkt;
    time_arr_opt(3, j) = toc;

end

for i = 1:num_methods
    t = time_arr(i,:);
    t(t == max(t) | t == min(t)) = [];
    times(i) = mean(t);

    t = time_arr_opt(i,:);
    t(t == max(t) | t == min(t)) = [];
    times_opt(i) = mean(t);
end

figure;
methods = {'GE', 'SOR', 'Schur-Newton'};
bar_data = [times, times_opt];
bar_colors = [1 0 1; 0 1 1]; % Magenta and cyan colors

% Create the bar plot
b = bar(bar_data, 'grouped');

% Set bar colors
for k = 1:numel(b)
    b(k).FaceColor = 'flat';
    b(k).CData = repmat(bar_colors(k,:), numel(methods), 1);
end

% Customize the graph
set(gca, 'XTickLabel', methods, 'XTick', 1:num_methods);
xlabel('Method');
ylabel('Execution Time (seconds)');
legend({'Unoptimized', 'Optimized'}, 'Location', 'northeast');
grid on;
grid minor;

%% Testing the symmetrized double-augmented IP method
