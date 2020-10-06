function [C, b_C] = ochs(dims, J, G, b_G)

% constants
kDimActualized      = dims.Actualized;
kDimUnActualized    = dims.UnActualized;

M = [zeros(kDimUnActualized, kDimActualized); eye(kDimActualized)];
U = null(J)';
U_bar = U*M;
[L, D] = ldl(U_bar'*U_bar);

TOL = 1e-8;
for i = 1:size(D,1)
    if abs(D(i,i)) < TOL
        D = D(:, 1:i-1);
        break;
    end
end
U_hat = (L*sqrt(D))';

C_bar = pinv(U_hat)';
C = [zeros(size(C_bar, 1), kDimUnActualized) C_bar];
% solve for b_C

JC = [J; C];
JCG = [JC; G];

if rank(JC) < rank(JCG)
    % infeasible goal
    C = [];
    b_C = [];
    return;
end

JG = [J; G];
b_JG = [zeros(size(J,1), 1); b_G];
% if rank([JG b_JG]) > rank(JG)
%     % infeasible problem
%     C = [];
%     b_C = [];
%     return;
% end

v_star = JG\b_JG;
b_C = C*v_star;
