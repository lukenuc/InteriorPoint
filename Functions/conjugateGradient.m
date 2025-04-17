
function [x, k] = conjugateGradient(A, b, M, x0, TOL, MAXIT)

x = x0(:); 
r = A*x0 - b; 
y = M.*r; 
p = -y; 

for k = 1:MAXIT
    if norm(r, inf) < TOL
        return
    else
        alpha = r'*y/(p'*A*p); 
        x = x + alpha*p; 
        rp1 = r + alpha*A*p;
        yp1 = M.*rp1; 
        beta = rp1'*yp1/(r'*y); 
        p = -yp1 + beta*p; 
        r = rp1; y = yp1; 
    end
end

end

