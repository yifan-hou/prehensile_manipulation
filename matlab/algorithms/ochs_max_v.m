% A Lambda <= b_A
function [action, time] = ochs_max_v(dims, J, G, b_G, Fg, A, b_A, J_All)
time.velocity = 0;
time.force = 0;
action.C = [];
tic
% constants
n_a = dims.Actualized;
n_u = dims.UnActualized;
n   = n_a + n_u;

M = [zeros(n_u, n_a); eye(n_a)];
U = null(J)';
rank_J = n - size(U,1);

if isempty(U)
    % no null space. The environment already fully constrained the problem
    return;
end

U_bar = U*M;

%% null space hfvc
U_bar_n = orth(U_bar')';
n_av_max  = size(U_bar_n,1);

JG = [J; G];
b_JG = [zeros(size(J,1), 1); b_G];
null_JG = null(JG);
rank_JG = n - size(null_JG,2);
if rank_JG - rank_J > n_av_max
    % infeasible problem: goal cannot be satisfied
    return;
end

n_av = n_av_max;
C_bar = U_bar_n;
C = [zeros(n_av, n_u) C_bar];

n_af  = n_a - n_av;
R_a = [null(C_bar)'; C_bar];
T = blkdiag(eye(n_u), R_a);

% solve for b_C
JC = [J; C];
JCG = [JC; G];

if rank(JC) < rank(JCG)
    % infeasible goal
    return;
end

% if rank([JG b_JG]) > rank(JG)
%     % infeasible problem
%     return;
% end

v_star = JG\b_JG;
b_C = C*v_star;

time_velocity = toc;


%% force part
tic;

nLambda = size(J_All,1);

Aeq_lp = [T*J_All' [zeros(n_u, n_a); eye(n_a)]];
beq_lp = -T*Fg;
A_lp = [A zeros(size(A,1), n_a)];
b_lp = b_A;

qp.Q = eye(nLambda + n_a);
qp.f = zeros(nLambda + n_a, 1);

options = optimoptions('quadprog', 'Display', 'none');
[x,~,EXITFLAG] = quadprog(qp.Q, qp.f, A_lp, b_lp, Aeq_lp, beq_lp, [], [], [],options);

if EXITFLAG == 1
    eta_af = x(nLambda + 1:nLambda + n_af);
    action.eta_af = eta_af;
end

action.n_av = n_av;
action.n_af = n_af;
action.R_a = R_a;
action.w_av = b_C;
action.C = C;
action.b_C = b_C;
action.qpflag = EXITFLAG;

time.velocity = time_velocity;
time.force = toc;