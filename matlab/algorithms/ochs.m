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
rank_J = n - size(U,1);

if isempty(U)
    % no null space. The environment already fully constrained the problem
    return;
end

U_bar = U*M;

%%
%% Minimize  cond(U'C')
%%
% U_hat = chol(U_bar'*U_bar + 1e-15*eye(n_a));
% TOL = 0.0001; % this tolerance has to be big.
% id = [];
% for i = 1:n_a
%     if (norm(U_hat(n_a-i + 1,:)) > TOL)
%         id = [id n_a-i+1];
%     end
% end
% U_hat = U_hat(id, :);
% U_hat = U_hat(id, :);
% n_av = size(U_hat,1);
% n_af = n_a - n_av;
% % C_bar = (U_hat\eye(n_av))';
% C_bar = lsqminnorm(U_hat, eye(n_av))';

%%
%% minimize cond(U'normalizeByRow(C'))
%%
% TOL = 0.0001; % this tolerance has to be big.
% id = [];
% row_norm = [];
% for i = 1:n_a
%     norm_i = norm(U_hat(n_a-i + 1,:));
%     if (norm_i > TOL)
%         id = [n_a-i+1 id];
%         row_norm = [norm_i row_norm];
%     end
% end
% over_dim = length(id) - size(U,1);
% if over_dim > 0
%     for i = 1:over_dim
%         [~, id_id] = min(row_norm);
%         id(id_id) = [];
%         row_norm(id_id) = [];
%     end
% end
% U_hat_null = null(U_hat);
% C_star = C_bar;
% if ~isempty(U_hat_null)
%     U_hat_null1 = U_hat_null(:,1);
%     cnorm = normByRow(C_bar);
%     [cn_max, cn_max_id] = max(cnorm);
%     for i = 1:n_av
%         if i == cn_max_id
%             continue;
%         end
%         k = sqrt(cn_max^2 - cnorm(i)^2);
%         C_star(i,:) = C_star(i,:) + k*U_hat_null1';
%     end
% end

%% null space hfvc
R_f = null(U_bar)';
R_v = null(R_f)';
n_av_max  = size(R_v,1);

JG = [J; G];
b_JG = [zeros(size(J,1), 1); b_G];
null_JG = null(JG);
rank_JG = n - size(null_JG,2);
if rank_JG - rankJ > n_av_max
    % infeasible problem: goal cannot be satisfied
    return;
end

% todo: what if rank_JG - rankJ = 0?
% may be redundant with other feasibility conditions

n_av = 0;
if rank_JG - rankJ == n_av_max
    % we have to use all free dof for velocity control
    n_av = n_av_max;
else
    % We don't need to use up free space for velocity control
    n_av = rank_JG - rankJ;
    null_JG
end

n_af  = n_a - n_av;


R_a = [R_f; R_v];
T = blkdiag(eye(n_u), R_a);
C = [zeros(n_av, n_u) R_v];

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