function solution = wrenchSpaceAnalysis_modeSelection(kFrictionE, kFrictionH, ...
            CP_W_e, CN_W_e, CP_H_h, CN_H_h, R_WH, p_WH, G, b_G, kForceMagnitude)

kCharacteristicLength = 0.05;
wrench_scaling_matrix = diag([1, 1, 1/kCharacteristicLength]);
TOL = 1e-7;

Ne = size(CP_W_e, 2);
Nh = size(CP_H_h, 2);

% add virtual z axis
assert(size(CP_W_e, 1) == 2);
CP_W_e = [CP_W_e; zeros(1, Ne)];
CN_W_e = [CN_W_e; zeros(1, Ne)];
CP_H_h = [CP_H_h; zeros(1, Nh)];
CN_H_h = [CN_H_h; zeros(1, Nh)];
R_WH = [R_WH [0; 0]; 0 0 1];
p_WH = [p_WH; 0];

R_HW = R_WH';
p_HW = -R_HW*p_WH;
adj_HW = SE22Adj(R_HW(1:2,1:2), p_HW(1:2));
adj_WH = SE22Adj(R_WH(1:2,1:2), p_WH(1:2));

[Jac_e, Jac_h] = getWholeJacobian(CP_W_e, CN_W_e, CP_H_h, CN_H_h, adj_WH, adj_HW);
[Jacf_e, Jacf_h] = getWholeJacobianFrictional(CP_W_e, CN_W_e, kFrictionE, ...
        CP_H_h, CN_H_h, kFrictionH, adj_WH, adj_HW);

fprintf("###############################################\n");
fprintf("##              Mode Enumeration             ##\n");
fprintf("###############################################\n");
% Enumerate contact modes
e_modes = contact_mode_enumeration(CP_W_e(1:2,:), CN_W_e(1:2,:), true);
h_modes = contact_mode_enumeration(CP_H_h(1:2,:), CN_H_h(1:2,:), true);

% get rid of all separation modes
e_modes(:, all(e_modes == 0, 1)) = [];
h_modes(:, all(h_modes == 0, 1)) = [];

% add all fixed mode if necessary
if sum(all(e_modes == 1, 1)) == 0
    e_modes = [e_modes ones(size(e_modes,1),1)];
end
if sum(all(h_modes == 1, 1)) == 0
    h_modes = [h_modes ones(size(h_modes,1),1)];
end
disp(['Enumerated modes: ' num2str(size(e_modes, 2)*size(h_modes, 2)) ' modes.']);

% EH cone intersection
%   * Compute safety margins
%   * get rid of infeasible cone
%

eh_modes = [];
margins = [];

eh_cones = cell(size(e_modes, 2)*size(h_modes, 2), 1);
cone_e = cell(size(e_modes, 2)*size(h_modes, 2), 1);
cone_h = cell(size(e_modes, 2)*size(h_modes, 2), 1);
Jacs = cell(size(e_modes, 2)*size(h_modes, 2), 1);
Jacus = cell(size(e_modes, 2)*size(h_modes, 2), 1);

eh_cone_feasible_mode_count = 0;
for i = 1:size(e_modes, 2)
    for j = 1:size(h_modes, 2)
        [Je_, Jh_] = getFrictionalJacobianFromContacts(e_modes(:, i), h_modes(:, j), Jacf_e, Jacf_h);

        % check force balance
        R = coneIntersection(Je_', Jh_');
        if isempty(R) || norm(R) == 0
            continue;
        end

        % compute velocity Jacobian
        [N, Nu] = getJacobianFromContacts(e_modes(:, i), h_modes(:, j), Jac_e, Jac_h);

        % compute cone stability margin
        % scaling
        Je_scaled_ = Je_*wrench_scaling_matrix';
        Jh_scaled_ = Jh_*wrench_scaling_matrix';
        R_scaled = wrench_scaling_matrix*R;

        margin_e = computeStabilityMargin(Je_scaled_', R_scaled);
        margin_h = computeStabilityMargin(Jh_scaled_', R_scaled);
        margin_ = min(margin_e, margin_h);

        % figure(1);clf(1);hold on;
        % printModes([e_modes(:, i); h_modes(:, j)]);
        % fprintf('Margin: %f\n', margin_);
        % drawWrench(Je_','g', true);
        % drawWrench(Jh_','b', true);

        eh_modes = [eh_modes [e_modes(:, i); h_modes(:, j)]];
        margins = [margins margin_];

        eh_cone_feasible_mode_count = eh_cone_feasible_mode_count + 1;
        eh_cones{eh_cone_feasible_mode_count} = R;
        Jacs{eh_cone_feasible_mode_count} = N;
        Jacus{eh_cone_feasible_mode_count} = Nu;
    end
end

disp(['Modes with none empty EH Cone: ' num2str(eh_cone_feasible_mode_count)]);

% trim
eh_cones = eh_cones(1:eh_cone_feasible_mode_count);
Jacs = Jacs(1:eh_cone_feasible_mode_count);
Jacus = Jacus(1:eh_cone_feasible_mode_count);

%%
%% Begin to check each cone
%%
fprintf("###############################################\n");
fprintf("##             Check Each Mode               ##\n");
fprintf("###############################################\n");

solutions = cell(eh_cone_feasible_mode_count, 1);
solutions_count = 0;

for m = 1:eh_cone_feasible_mode_count
    disp('***********');
    fprintf("Mode %d of %d\n", m, eh_cone_feasible_mode_count);
    eh_mode_goal = eh_modes(:, m);
    printModes(eh_mode_goal);

    fprintf("=======================\n");
    fprintf("=== Hybrid Servoing ===\n");
    fprintf("=======================\n");
    %   velocity vector: v = [vo, vh]
    % goal velocity
    dims.Actualized = 3;
    dims.UnActualized = 3;
    dims.SlidingFriction = 0;

    N_all = Jacs{m};

    [n_av, n_af, R_a, w_av] = solvehfvc(dims, N_all, G, ...
        b_G, [], [], [], [], [], 3, false);

    if isempty(R_a) || any(isnan(w_av))
        disp('Hybrid servoing has no solution.')
        continue;
    end
    fprintf('The force controlled dimension is %d\n', n_af);
    fprintf('The velocity controlled dimension is %d\n', n_av);
    assert(n_af < 3);
    assert(n_af > 0);

    % make sure all velocity commands >= 0
    R_id_flip = find(w_av < 0);
    R_a(n_af + R_id_flip, :) = - R_a(n_af + R_id_flip, :);
    w_av(R_id_flip) = - w_av(R_id_flip);
    R_a_inv = R_a^-1;

    Cv = [zeros(n_av, 3), R_a(n_af + 1: end, :)];
    b_C = w_av;

    fprintf("=======================\n");
    fprintf("=== Mode Filtering ===\n");
    fprintf("=======================\n");

    % How to filter out modes:
    % 1. If NC degenerates, mark this mode as incompatible;
    % 2. If nominal velocity under NC exists, and it cause the contact point to
    %    slide in a different direction, remove this mode.
    feasible_mode_id = false(eh_cone_feasible_mode_count, 1);

    ehw_cones_m = cell(length(eh_cone_feasible_mode_count), 1);

    % cone of possible actuation forces
    %   force controlled direction can have both positive and negative forces
    W_action = [R_a_inv(:, 1:n_af), -R_a_inv];

    TOL = 1e-7;
    feasible_mode_count = 0;
    flag_crashing = false;
    flag_goal_infeasible = false;
    flag_goal_impossible = false;
    flag_goal_w_impossible = false;

    for n = 1:eh_cone_feasible_mode_count
        % filter out modes using velocity command
        N = Jacs{n};
        Nu = Jacus{n};

        N = rref(N);
        rank_N = rank(N);
        N = N(1:rank_N, :);

        % compute possible sliding directions
        % these computations are based on 'Criteria for Maintaining Desired Contacts for Quasi-Static Systems'
        Lambda_bar = [Cv; N];
        b_Lambda_bar = [b_C; zeros(size(N,1), 1)];

        compatible = false;
        if rank([N;Cv], TOL) - rank_N > 0
            v_star = linsolve(Lambda_bar, b_Lambda_bar);
            if any(Nu*v_star < -TOL)
                % this mode can not exist
                % V-Impossible
                if n == m
                    flag_goal_impossible = true;
                    break;
                else
                    continue;
                end
            end
            compatible = true;
        else
            if n == m
                flag_goal_infeasible = true;
                break;
            end
        end

        % check for feasible actuation
        R_ = coneIntersection(eh_cones{n}, W_action);
        if isempty(R_) || norm(R_) == 0
            if n == m
                flag_goal_w_impossible = true;
                break;
            else
                continue;
            end
        end

        % Crashing check:
        % Check if the infeasible mode contains v-controlled direction
        if ~compatible
            intersection = coneIntersection(R_, W_action(:, end-n_av+1:end));
            if ~isempty(intersection) && norm(intersection) > TOL
                flag_crashing = true;
                break;
            else
                continue;
            end
        end

        % figure(1);clf(1);hold on;
        % printModes([e_modes(:, i); h_modes(:, j)]);
        % fprintf('Margin: %f\n', margin_);
        % drawWrench(Je_','g', true);
        % drawWrench(Jh_','b', true);
        % drawWrench(W_action,'k', true);

        feasible_mode_count = feasible_mode_count + 1;
        feasible_mode_id(n) = true;
        ehw_cones_m{feasible_mode_count} = R_;
    end

    if flag_crashing
        disp('Crashing. Discard this mode.');
        continue;
    end
    if flag_goal_infeasible
        disp('Goal mode violates velocity equality constraints. Discard this mode.');
        continue;
    end
    if flag_goal_impossible
        disp('Goal mode violates velocity inequality constraints. Discard this mode.');
        continue;
    end

    if flag_goal_w_impossible
        disp('Goal mode force can not be maintained by velocity. Discard this mode.');
        continue;
    end

    disp(['Remaining feasible modes: ' num2str(feasible_mode_count)]);

    % trim
    feasible_mode_id = find(feasible_mode_id);

    goal_id = find(feasible_mode_id == m);
    others_id = find(feasible_mode_id ~= m);

    eh_cones_goal_m = ehw_cones_m{goal_id};
    eh_cones_other_feasible_m = ehw_cones_m(others_id);

    fprintf("=======================\n");
    fprintf("===   Evaluation    ===\n");
    fprintf("=======================\n");

    % action selection
    flag_force_region_feasible = true;

    % project to force controlled plane/line
    force_basis = R_a_inv(:, 1:n_af);
    projection_goal = force_basis'*eh_cones_goal_m;
    projection_goal_remains = projection_goal;

    force_action = mean(normc(projection_goal), 2);
    shape_margin = [];
    if ~isempty(eh_cones_other_feasible_m)
        if n_af == 1
            shape_margin = norm(wrench_scaling_matrix*force_action);
            for i = 1:length(eh_cones_other_feasible_m)
                projection_other_feasible = force_basis'*eh_cones_other_feasible_m{i};
                if any(projection_other_feasible'*force_action > 0)
                    flag_force_region_feasible = false;
                    break;
                end
            end
        elseif n_af == 2
            for i = 1:length(eh_cones_other_feasible_m)
                projection_other_feasible = force_basis'*eh_cones_other_feasible_m{i};
                [~, projection_goal_remains] = coneIntersection2D( ...
                        projection_goal_remains, projection_other_feasible);
                if isempty(projection_goal_remains)
                    flag_force_region_feasible = false;
                    break;
                end
            end
            if flag_force_region_feasible
                % pick the action that is far away from other regions.
                if length(eh_cones_other_feasible_m) > 1
                    force_action = mean(normc(projection_goal_remains), 2);
                else
                    ang1 = angBTVec(projection_other_feasible, projection_goal_remains(:, 1));
                    ang2 = angBTVec(projection_other_feasible, projection_goal_remains(:, 2));
                    if ang1 > ang2
                        force_action = normc(projection_goal_remains(:, 1));
                    else
                        force_action = normc(projection_goal_remains(:, 2));
                    end
                end
                projection_goal_remains_scaled = wrench_scaling_matrix*force_basis*projection_goal_remains;
                shape_margin = kForceMagnitude*angBTVec(projection_goal_remains_scaled(:, 1), projection_goal_remains_scaled(:, 2))/2;
            end
        end
    end

    if flag_force_region_feasible
        disp('Force command:');
        disp(force_action);
        solution.eh_mode = eh_mode_goal;
        solution.n_af = n_af;
        solution.n_av = n_av;
        solution.R_a = R_a;
        solution.w_av = w_av;
        solution.eta_af = -kForceMagnitude*force_action; % the minus sign comes from force balance
        solution.margin = min(shape_margin, margins(m));

        solutions_count = solutions_count + 1;
        solutions{solutions_count} = solution;
    else
        disp('No distinguishable force command.');
        continue;
    end
end
fprintf("###############################################\n");
fprintf("##                  Results                  ##\n");
fprintf("###############################################\n");
fprintf("Total number of solutions found: %d\n", solutions_count);
solution = [];
if solutions_count > 0
    solutions = solutions(1:solutions_count);
    margins = zeros(1, solutions_count);
    for i = 1:solutions_count
        printModes(solutions{i}.eh_mode);
        fprintf("margin: %f\n", solutions{i}.margin);
        margins(i) = solutions{i}.margin;
    end

    [~, best_solution_id] = max(margins);
    solution = solutions{best_solution_id};
end


return;
