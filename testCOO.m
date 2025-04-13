
clear; clc

% A = [1 4 2; 4 -1 3; 6 3 8]; 
% b = [1; 2; -1]; 
n = 100; rho = 0.1; 
A = eye(n) + full(sprand(n, n, rho)); 
b = rand(n,1); 

[row, col, val] = find(A); 
tic
[urow, ucol, uval, c, exitflag] = GaussElimCOO(row, col, val, b, n); 
toc
U = sparse(urow, ucol, uval); 
U = full(U); 
x_sol = backSubs(U,c)

% tic
% [U,c] = GaussElimPivoting(A,b); 
% backSubs(U,c);
% toc