function [Minv] = schurInverseNewton(M, p, q)

persistent MsDinv_prev;

% this all needs to be updated for the swapped rows of the gradient
A = M(1:p,1:p); B = M(1:p,p+1:end);
C = M(p+1:end,1:p); D = M(p+1:end,p+1:end);

Z = D(q+1:end, 1:q); S = D(q+1:end, q+1:end);
diag_sinv = 1 ./ diag(S);
diag_z = diag(Z);

Dinv = [eye(q) zeros(q); diag(-diag_sinv.*diag_z) diag(diag_sinv)];
MsD = A - B*Dinv*C;
% MsDinv = inv(MsD);

if(~isempty(MsDinv_prev))
    MsDinv = mtx_inv_newton(MsD, 1.e-6, 100, MsDinv_prev);
else
    MsDinv = mtx_inv_newton(MsD, 1.e-6, 100, MsD'/(norm(MsD,1)*norm(MsD,inf)));
end

T = MsDinv*B*Dinv; U = Dinv*C;
Minv = [MsDinv -T; -U*MsDinv Dinv+U*T];
MsDinv_prev = MsDinv; 

end

