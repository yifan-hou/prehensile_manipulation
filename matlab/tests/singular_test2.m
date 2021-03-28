clear;clc;
p_We = [0;0];
n_We = [0;1];

p_Hh = [0;
        0];

n_Hh = [0;
        -1];

R_WH = eye(2);
p_WH = [0; 0.2];

G = [1 0 0 0 0 0];
b_G = 0.1;

Fg = zeros(6,1);
Fg(2) = -10;
    
num_seeds = 5;

R_HW = R_WH';
p_HW = -R_HW*p_WH;
adj_HW = SE22Adj(R_HW, p_HW);
adj_WH = SE22Adj(R_WH, p_WH);

kFrictionE = 0.3;
kFrictionH = 0.8;
emodes = [1];
hmodes = [2];
[N_e, T_e, N_h, T_h, eCone, eTCone, hCone, hTCone] = getWholeJacobian(p_We, n_We, ...
    p_Hh, n_Hh, adj_WH, adj_HW, 1, kFrictionE, kFrictionH);
[J, ~, normal_ids] = getJacobianFromContacts(emodes, hmodes, N_e, N_h, T_e, T_h);

nLambda = size(J,1);
A = eye(nLambda);
A = -A(normal_ids, :);
b_A = -5*ones(size(A,1),1);

%
dims.Actualized = 3;
dims.UnActualized = 3;

[action, time] = ochs(dims, J, G, b_G, Fg, A, b_A, J);
C1 = action.C;
[action, time] = hybrid_servoing(dims, J, G, b_G, Fg, A, b_A, num_seeds, J);
C2 = action.C;
M = blkdiag(zeros(dims.UnActualized), eye(dims.Actualized));
U = null(J)';
U_bar = U*M;
C3 = null(null(U_bar)')';

Jreg = rref(J);
Jreg = Jreg(1:rank(J), :);
score1 = cond(normalizeByRow([Jreg; C1]));
score2 = cond(normalizeByRow([Jreg; C2]));
score3 = cond(normalizeByRow([Jreg; C3]));

disp([score1 score2 score3]);
    