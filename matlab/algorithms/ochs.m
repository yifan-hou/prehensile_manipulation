% A Lambda <= b_A
function [C, b_C, time] = ochs(dims, J, G, b_G, Fg, A, b_A)
time.velocity = 0;
time.force = 0;
tic
% constants
n_a = dims.Actualized;
n_u = dims.UnActualized;
n   = n_a + n_u;

M = [zeros(n_u, n_a); eye(n_a)];
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

time_velocity = toc;


%% force part
tic;

M = [T*J' [zeros(n - n_av, n_av); eye(n_av)]];
n_free = size(M,2);
n_dual_free = size(M,1);

Aeq_lp = [2*eye(n_free), M', zeros(n_free, n_af);
          M, zeros(n_dual_free), [zeros(n_u, n_af); eye(n_af); zeros(n_av, n_af)]];
beq_lp = [zeros(n_free, 1); -T*Fg];
A_lp = [A zeros(size(A,1), size(Aeq_lp,2) - size(A,2))];
b_lp = b_A;


qp.Q = diag([zeros(1, n_free + n_dual_free), ones(1, n_af)]);
qp.f = zeros(n_free + n_dual_free + n_af, 1);

if (n_af == 0)
    qp.Q = eye(n_free + n_dual_free + n_af);
end

options = optimoptions('quadprog', 'Display', 'none');
x = quadprog(qp.Q, qp.f, A_lp, b_lp, Aeq_lp, beq_lp, [], [], [],options);


eta_af = x(n_free + n_dual_free + 1:end);

time.velocity = time_velocity;
time.force = toc;