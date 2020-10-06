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
%   TODO: the use of dims.slidingfriction is weird
%         do rref on N when solving for v_star
%   TODO in c++:
%         handle n_av_min=0 return in c++ code

function [C, b_C] = hybrid_servoing(dims, J, G, b_G, num_seeds)
kNumSeeds = num_seeds;

% constants
kDimActualized      = dims.Actualized;
kDimUnActualized    = dims.UnActualized;

kDimLambda       = size(J, 1);
kDimGeneralized  = kDimActualized + kDimUnActualized;

JG = [J; G];

rank_J = rank(J);
rank_JG = rank(JG);

n_av_min = rank_JG - rank_J;

n_av = n_av_min;
if n_av_min == 0
    % infeasible goal
    C = [];
    b_C = [];
    return;
end

% b_JG = [zeros(size(N, 1), 1); b_G];
JG_nullspace_basis = null(JG);
basis_c = null([JG_nullspace_basis';
        eye(kDimUnActualized), zeros(kDimUnActualized,kDimActualized)]);

if (rank_J + kDimActualized < kDimGeneralized)
    % not prehensile, infeasible
    C = [];
    b_C = [];
    return;
end

n_c = rank_JG - kDimUnActualized;
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
    C = [];
    b_C = [];
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

R_a = [null(c_best(:, kDimUnActualized+1:end))';
        c_best(:, kDimUnActualized+1:end)];
T = blkdiag(eye(kDimUnActualized), R_a);

w_av = c_best*v_star;

C = T(end - n_av + 1:end, :);
b_C = w_av;