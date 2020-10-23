%% Solve the hybrid servoing problem. The algorithm solves for a hybrid force-velocity control command described by
% dimensionality @p n_av, @p n_af, direction @p R_a and magnitudes @p w_av, @p eta_af.
%
% The constraints are natural holonomic constraints: @p N_all * v = 0, goal constraints @p G*v = @p b_G, static force balance, and contact mode guard conditions.
%
% The system has kDOF = dims.Actualized +
% dims.UnActualized dimensions. This means the generalized force f and
% generalized velocity v are kDOF dimensional.
% The first dims.UnActualized dimensions correspond to the free DOFs in the
% system (such as free objects), while the last dims.Actualized DOFs correspond
% to the controllable DOFs (such as robot joints.)
%
%
% N_all: Linear constraints on generalized velocity. N_all v = 0
% G, b_G: Goal description, affine constraints on generalized velocity
%         G*v = b_G
% F: External force vector. Same size as generalized force
% Aeq, beq: Guard condition, Aeq*v = beq
% A, b_A: Guard condition, A*v <= b_A
% dims:
%     dims.Actualized: number of actualized dimensions
%     dims.UnActualized: number of unactualized dimensions
%     dims.SlidingFriction: number of sliding friction dimensions
%
%   optional:
%       num_seeds: number of random initializations to try when solving for
%               the velocity control
% Outputs
%   n_av: number, dimensionality of velocity controlled actions
%   n_af: number, dimensionality of force controlled actions
%   R_a: (n_av+n_af)x(n_av+n_af) matrix, the transformation that describes the
%           direction of velocity&force control actions. f first, v latter.
%   w_av: n_av x 1 vector, magnitudes of velocity controls
%   eta_af: n_af x 1 vector,  magnitudes of force controls
%
%   TODO: use J_All in force computation
function [action, time] = hybrid_servoing(dims, J, G, b_G, Fg, A, b_A, num_seeds, J_All)
time.velocity = 0;
time.force = 0;
action.C = [];
tic;
kNumSeeds = num_seeds;

% constants
n_a = dims.Actualized;
n_u = dims.UnActualized;

kDimLambda = size(J, 1);
n = n_a + n_u;

JG = [J; G];

rank_J = rank(J);
rank_JG = rank(JG);

n_av_min = rank_JG - rank_J;

n_av = n_av_min;
n_af = n_a - n_av;
if (n_af == 0)
    % no need for HFVC; all velocity control. Can not satisfy guard
    % condition
    return;
end
if n_av_min == 0
    % infeasible goal
    return;
end

% b_JG = [zeros(size(N, 1), 1); b_G];
JG_nullspace_basis = null(JG);
basis_c = null([JG_nullspace_basis';
        eye(n_u), zeros(n_u,n_a)]);


if (rank_J + n_a < n)
    % not prehensile, infeasible
    return;
end

n_c = rank_JG - n_u;
if (size(basis_c, 2) > n_c)
    % has uncontrollable DOF, infeasible
    return;
end
b_JG = [zeros(kDimLambda, 1); b_G];
% if rank([JG b_JG]) > rank_JG
%     % Goal is infeasible
%     C = [];
%     b_C = [];
%     return;
% end
v_star = JG\b_JG;

% Projected gradient descent
NIter             = 50;
J_nullspace_basis = null(J);
BB                = basis_c'*basis_c;
JJ                = J_nullspace_basis*(J_nullspace_basis');

if norm((basis_c')*JJ*basis_c) < 1e-9
    % can not make C out of N, no solution
    return;
end
% added a term:
%  min kCoefNu*c_i'*v_star*N_u*c_i

% VstarN_u = normc(v_star)*sum(normr(N_u));
% kCoefNu = 0;
cost_all = zeros(1, kNumSeeds);
k_all = rand([n_c, n_av, kNumSeeds]);
for seed = 1:kNumSeeds
    k  = k_all(:,:, seed);
    kn = normByCol(basis_c*k);
    k  = bsxfun(@rdivide, k, kn);

    for iter = 1:NIter
        % compute gradient0
        g = zeros(n_c, n_av);
        costs = 0;
        for i = 1:n_av
            ki = k(:,i);
            for j = 1:n_av
                if i == j
                    continue;
                end
                kj = k(:,j);
                costs = costs + (ki'*BB*kj)^2;
                g(:, i) = g(:, i) + 2*(ki'*BB*kj)*BB*kj;
            end
            g(:, i) = g(:, i) - 2*(basis_c')*JJ*basis_c*ki;
            % g(:, i) = g(:, i) + kCoefNu*2*(basis_c')*VstarN_u*basis_c*ki;
            costs   = costs - ki'*(basis_c')*JJ*basis_c*ki;
            % costs   = costs + kCoefNu*ki'*(basis_c')*VstarN_u*basis_c*ki;
        end
        % disp(['cost: ' num2str(costs)]);
        % descent
        delta = 0.1;
        k     = k - delta*g;
        % project
        kn = normByCol(basis_c*k);
        k  = bsxfun(@rdivide, k, kn);
    end
    cost_all(seed) = costs;
    k_all(:,:,seed) = k;
    % disp(['cost: ' num2str(costs)]);
end

[~, best_id] = min(cost_all);
k_best = k_all(:,:,best_id);
c_best = (basis_c*k_best)';

R_a = [null(c_best(:, n_u+1:end))';
        c_best(:, n_u+1:end)];
T = blkdiag(eye(n_u), R_a);

w_av = c_best*v_star;

C = T(end - n_av + 1:end, :);
b_C = w_av;

time_velocity = toc;




% Force
tic;
% unactuated dimensions
H = [eye(n_u), zeros(n_u, n_a)];
% Newton's laws
T_inv = T^-1;
M_newton = [zeros(n_u, kDimLambda) H*T_inv; ...
            T*J' eye(n)];
b_newton = [zeros(size(H,1), 1); -T*Fg];

M_free = M_newton(:, [1:kDimLambda+n_u, kDimLambda+n_u+n_af+1:end]);
M_eta_af = M_newton(:, [kDimLambda+n_u+1:kDimLambda+n_u+n_af]);

% prepare the QP
%   variables: [free_force, dual_free_force, eta_af]
n_free = kDimLambda + n_u + n_av;
n_dual_free = size(M_newton, 1);
% 0.5 x'Qx + f'x
qp.Q = diag([zeros(1, n_free + n_dual_free), ones(1, n_af)]);
qp.f = zeros(n_free + n_dual_free + n_af, 1);

% Aeq = beq
qp.Aeq = [2*eye(n_free), M_free', zeros(n_free, n_af);
          M_free, zeros(size(M_free, 1)), M_eta_af];
qp.beq = [zeros(n_free, 1); b_newton];

% Ax<b
qp.A = [A zeros(size(A,1), size(qp.Aeq,2) - size(A,2))];
qp.b = b_A;

options = optimoptions('quadprog', 'Display', 'none');
x = quadprog(qp.Q, qp.f, qp.A, qp.b, qp.Aeq, qp.beq, [], [], [],options);

eta_af = x(n_free + n_dual_free + 1:end);

action.n_av = n_av;
action.n_af = n_af;
action.R_a = R_a;
action.eta_af = eta_af;
action.w_av = w_av;
action.C = C;
action.b_C = b_C;

time.velocity = time_velocity;
time.force = toc;


