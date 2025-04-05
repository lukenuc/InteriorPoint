function [Minv] = schurInverse(M, p, q)

A = M(1:p,1:p); B = M(1:p,p+1:end);
C = M(p+1:end,1:p); D = M(p+1:end,p+1:end);

Z = D(1:q, 1:q); S = D(1:q, q+1:end); 
diag_sinv = 1 ./ diag(S); 
diag_z = diag(Z); 

Dinv = [zeros(q) eye(q); diag(diag_sinv) diag(-diag_sinv.*diag_z)]; 
MsD = A - B*Dinv*C; 
MsDinv = inv(MsD);

T = MsDinv*B*Dinv; U = Dinv*C;
Minv = [MsDinv -T; -U*MsDinv Dinv+U*T];
end

