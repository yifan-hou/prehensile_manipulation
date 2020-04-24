% Jac includes contact screws for all friction cone edges.
% Nv >= 0
% N = [Je, 0; -Jh Jh]
function [Je, Jh, CNF_e, CNF_h] = getWholeJacobianFrictional(CP_W_e, CN_W_e, ...
        mu_e, CP_H_h, CN_H_h, mu_h, adj_WH, adj_HW)

% friction cone
[CPF_e, CNF_e] = frictionCone2D(CP_W_e, CN_W_e, mu_e);
[CPF_h, CNF_h] = frictionCone2D(CP_H_h, CN_H_h, mu_h);

% contact screws
W_h_all = contactScrew(CPF_h, CNF_h);
W_e_all = contactScrew(CPF_e, CNF_e);

% scale based on friction coefficient
W_h_all = (1 + mu_h^2) * W_h_all;
W_e_all = (1 + mu_e^2) * W_e_all;

% the linear constraint on velocity from contact screws
Je = v3t2(W_e_all)'*adj_WH;
Jh = -v3t2(W_h_all)'*adj_HW*adj_WH;