function [Minv] = schurWoodbury(M, p, q)

persistent A; 
persistent Ainv; 
persistent B; 
persistent C; 
persistent MsA; 

D = M(p+1:end,p+1:end);
if isempty(A)
    A = M(1:p,1:p); 
    Ainv = inv(A); 
    B = M(1:p,p+1:end);
    C = M(p+1:end,1:p);
    prev_diag_s = zeros(q,1); 
    prev_diag_z = zeros(q,1); 
    MsAinv = inv(D - C*Ainv*B); 
end

Z = D(q+1:end, 1:q); S = D(q+1:end, q+1:end);
diag_s = diag(S);
diag_z = diag(Z);
U = [zeros(q); eye(q)]; % figure out dimensions here
V = [diag(diag_z-prev_diag_z) diag(diag_s-prev_diag_s)];

T = V*MsAinv; 
MsAinv = MsAinv - MsAinv*U*(inv(eye(q) + T*U))*T; 

T = MsAinv*C*Ainv; U = Ainv*B;
Minv = [Ainv + U*T -U*MsAinv; -T MsAinv];

%% method using M/A

% persistent A; 
% persistent Ainv; 
% persistent B; 
% persistent C; 
% persistent prev_diag_s; 
% persistent prev_diag_z; 
% persistent MsAinv; 
% 
% D = M(p+1:end,p+1:end);
% if isempty(A)
%     A = M(1:p,1:p); 
%     Ainv = inv(A); 
%     B = M(1:p,p+1:end);
%     C = M(p+1:end,1:p);
%     prev_diag_s = zeros(q,1); 
%     prev_diag_z = zeros(q,1); 
%     MsAinv = inv(D - C*Ainv*B); 
% end
% 
% Z = D(q+1:end, 1:q); S = D(q+1:end, q+1:end);
% diag_s = diag(S);
% diag_z = diag(Z);
% U = [zeros(q); eye(q)]; % figure out dimensions here
% V = [diag(diag_z-prev_diag_z) diag(diag_s-prev_diag_s)];
% 
% T = V*MsAinv; 
% MsAinv = MsAinv - MsAinv*U*(inv(eye(q) + T*U))*T; 
% 
% T = MsAinv*C*Ainv; U = Ainv*B;
% Minv = [Ainv + U*T -U*MsAinv; -T MsAinv];

end

