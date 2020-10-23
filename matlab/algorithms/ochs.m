% A Lambda <= b_A
function [action, time] = ochs(dims, J, G, b_G, Fg, A, b_A, J_All)
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

if isempty(U)
    % no null space. The environment already fully constrained the problem
    return;
end

U_bar = U*M;
U_hat = chol(U_bar'*U_bar + 1e-15*eye(n_a));

TOL = 0.0001; % this tolerance has to be big.
% for i = 1:n_a
%     if (norm(U_hat(n_a-i + 1,:)) > TOL)
%         U_hat = U_hat(1:n_a-i+1, :);
%         break;
%     end
% end
id = [];
for i = 1:n_a
    if (norm(U_hat(n_a-i + 1,:)) > TOL)
        id = [id n_a-i+1];
    end
end
U_hat = U_hat(id, :);


n_av = size(U_hat,1);
n_af = n_a - n_av;

C_bar = (U_hat\eye(n_av))';
R_a = [null(C_bar)';
        C_bar];
T = blkdiag(eye(n_u), R_a);
C = [zeros(n_av, n_u) C_bar];
% solve for b_C

JC = [J; C];
JCG = [JC; G];

if rank(JC) < rank(JCG)
    % infeasible goal
    return;
end

JG = [J; G];
b_JG = [zeros(size(J,1), 1); b_G];
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
x = quadprog(qp.Q, qp.f, A_lp, b_lp, Aeq_lp, beq_lp, [], [], [],options);

eta_af = x(nLambda + 1:nLambda + n_af);

action.n_av = n_av;
action.n_af = n_af;
action.R_a = R_a;
action.eta_af = eta_af;
action.w_av = b_C;
action.C = C;
action.b_C = b_C;

time.velocity = time_velocity;
time.force = toc;