% Jac includes contact screws for all friction cone edges.
% Nv >= 0
% N = [Jac_e, 0; -Jac_h Jac_h]
function [Jac_e, Jac_h] = getWholeJacobianFrictional(CP_e, CN_e, mu_e, CP_O_h, CN_O_h, mu_h, Adj_OW)

% friction cone
[CPF_e, CNF_e] = frictionCone2D(CP_e, CN_e, mu_e);
[CPF_h, CNF_h] = frictionCone2D(CP_O_h, -CN_O_h, mu_h); % note the minus sign

% contact screws
W_h_all = contactScrew(CPF_h, CNF_h);
W_e_all = contactScrew(CPF_e, CNF_e);

% the linear constraint on velocity from contact screws
Jac_e = v3t2(W_e_all)';
Jac_h = v3t2(W_h_all)'*Adj_OW;