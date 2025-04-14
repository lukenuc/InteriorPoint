% This function uses SOR's method to solve a linear system with
% matrix A and rhs b.
% Can call this function as [x,k,exitflag] = SOR(A,b,x0,omega,TOL,MAXIT)
% Inputs:
%   - A:     square nxn matrix
%   - b:     nx1 rhs vector
%   - x0:    initial guess
%   - TOL:   tolerance
%   - MAXIT: maximum number of iterations
% Outputs:
%   - x: solution (assuming the method converged)
%   - k: number of iterations
%   - exitflag: = -1 the method did not converge,
%               =  0 the method converged.
%
function [x,k,exitflag] = SOR_COO(arow,acol,aval,b,x0,omega,TOL,MAXIT)

x = x0(:);
n = length(b);
l = length(aval); 

for k = 1 : MAXIT
    y = x; % Save previous iteration
    cnt = 1; 
    for i = 1 : n
        sum1 = 0; sum2 = 0; aii = 0; 
        if i == 141
            debug = 1; 
        end
        while (cnt <= l & arow(cnt) == i)
            col = acol(cnt); a = aval(cnt);
            if col < i 
                sum1 = sum1 + a*x(col);
            elseif col > i                 
                sum2 = sum2 + a*y(col); 
            else
                aii = a; 
            end
            cnt = cnt + 1; 
        end
        x(i) = (omega*b(i) - omega*sum1 - omega*sum2 + (1-omega)*aii*y(i))/aii;
    end

    if norm(y-x,inf) < TOL
        exitflag = 0;
        return
    end

end

exitflag = -1;