
% =========================================================================
% FUNCTION NAME: GaussElimPivoting.m
% AUTHOR:        Luke Nuculaj
% CLASS:         APM 5334 - Applied Numerical Methods
% DESCRIPTION:  This function performs Gaussian elimination on a matrix A and rhs b with
% partial pivoting for numerical stability.
%
% INPUTS:
%   - A: square nxn matrix
%   - b: nx1 rhs vector
%
% OUTPUTS:
%   - A: upper triangular form of A
%   - b: modified rhs
%
% EXIT FLAGS:
%   - exitflag  = -1 means A is not square,
%               = -2 means b is not a column vector or size not
%                    compatible with A,
%               = 0 means the algorithm was successful.
% =========================================================================

function [A,b,exitflag] = GaussElim_IP(A,b,p,q)

n = size(A,1); % number of rows
m = size(A,2); % number of columns
tol = 1.e-12;

% Check that A is square
if n ~= m
    exitflag = -1;
    return
end

% Check b is column vector and compatible with A
if n ~= size(b,1) || size(b,2) ~= 1
    exitflag = -2;
    return
end

for j = 1 : n-1
    % spy(abs(A) > tol)
    if (j <= p) || (j >= p+q+1)
        for i = j+1 : n
            m = A(i,j)/A(j,j);
            for k = j : n
                A(i,k) = A(i,k) - m*A(j,k);
            end
            b(i) = b(i) - m*b(j);
        end
    else
        m = 1/A(j,j);
        A(j+q,j) = 0;
        A(j+q,j+q) = A(j+q,j+q) - m*A(j,j+q);
        b(j+q) = b(j+q) - m*b(j);
    end
end

exitflag = 0;