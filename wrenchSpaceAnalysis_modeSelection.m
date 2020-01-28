function solution = wrenchSpaceAnalysis_modeSelection(...
        kFrictionE, kFrictionH, CP_W_e, CN_W_e,...
        CP_H_h, CN_H_h, R_WH, p_WH, G, b_G, kForceMagnitude)
solution = [];
TOL = 1e-7;

% scaling for generalized velocity
% V = gvscale * V_scaled
kCharacteristicLength = 0.15;
vscale = diag([1 1 1/kCharacteristicLength]);
vscale_inv = diag([1 1 kCharacteristicLength]);
gvscale = diag([1 1 1/kCharacteristicLength 1 1 1/kCharacteristicLength]);
% N_ * V_scaled = 0, N*V = 0,
% -> N_ = N*gvscale

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

% scale
Jac_e = Jac_e * gvscale;
Jac_h = Jac_h * gvscale;
Jacf_e = Jacf_e * vscale;
Jacf_h = Jacf_h * vscale;

R_all_f = coneIntersection(Jacf_e', Jacf_h');

fprintf("###############################################\n");
fprintf("##              Mode Enumeration             ##\n");
fprintf("###############################################\n");
[e_modes, h_modes] = partialGraspModeEnumeration(CP_W_e, CN_W_e, CP_H_h, CN_H_h);

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
Jacues = cell(size(e_modes, 2)*size(h_modes, 2), 1);

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
        [N, Nu, Nue] = getJacobianFromContacts(e_modes(:, i), h_modes(:, j), Jac_e, Jac_h);

        % compute cone stability margin
        if sameConeCheck(Je_', R)
            margin_ = coneStabilityMargin(Jh_', R);
        elseif sameConeCheck(Jh_', R)
            margin_ = coneStabilityMargin(Je_', R);
        else
            margin_e = coneStabilityMargin(Je_', R);
            margin_h = coneStabilityMargin(Jh_', R);
            margin_ = min(margin_e, margin_h);
        end

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
        Jacues{eh_cone_feasible_mode_count} = Nue;
    end
end

disp(['Modes with none empty EH Cone: ' num2str(eh_cone_feasible_mode_count)]);

% trim
eh_cones = eh_cones(1:eh_cone_feasible_mode_count);
Jacs = Jacs(1:eh_cone_feasible_mode_count);
Jacus = Jacus(1:eh_cone_feasible_mode_count);
Jacues = Jacues(1:eh_cone_feasible_mode_count);

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

    fprintf("====================================\n");
    fprintf("= Hybrid Servoing & Crashing Check =\n");
    fprintf("====================================\n");
    N_all = Jacs{m};
    N_ue = Jacues{m};
    [n_av, n_af, R_a, R_a_inv, w_av,Cv, b_C] = hybridServoing(N_all, N_ue, G, b_G);
    if isempty(n_av)
        continue;
    end
    W_action = [R_a_inv(:, 1:n_af), -R_a_inv];
    intersection = coneIntersection(R_all_f, W_action(:, end-n_av+1:end));
    % Crashing check
    if ~isempty(intersection) && norm(intersection) > TOL
        disp('Crashing. Discard this mode.');
        continue;
    end

    fprintf("=======================\n");
    fprintf("=== Mode Filtering ===\n");
    fprintf("=======================\n");

    % How to filter out modes:
    % 1. If NC degenerates, mark this mode as incompatible;
    % 2. If nominal velocity under NC exists, and it cause the contact point to
    %    slide in a different direction, remove this mode.
    feasible_mode_id = false(eh_cone_feasible_mode_count, 1);

    % cone of possible actuation forces
    %   force controlled direction can have both positive and negative forces

    TOL = 1e-7;
    feasible_mode_count = 0;
    flag_goal_infeasible = false;
    flag_goal_impossible = false;


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
        if rank([Lambda_bar b_Lambda_bar], TOL) == rank(Lambda_bar, TOL)
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

        % figure(1);clf(1);hold on;
        % drawWrench(eh_cones{n},'g', true);
        % drawWrench(W_action,'k', true);
        if compatible
            feasible_mode_count = feasible_mode_count + 1;
            feasible_mode_id(n) = true;
        end
    end

    if flag_goal_infeasible
        disp('Goal mode violates velocity equality constraints. Discard this mode.');
        continue;
    end
    if flag_goal_impossible
        disp('Goal mode violates velocity inequality constraints. Discard this mode.');
        continue;
    end
    disp(['Remaining feasible modes: ' num2str(feasible_mode_count)]);

    eh_cones_goal_m = eh_cones{m};
    feasible_mode_id(m) = 0;
    eh_cones_other_feasible_m = eh_cones(feasible_mode_id);

    fprintf("===============================\n");
    fprintf("===  Compute Force Action   ===\n");
    fprintf("===============================\n");

    % action selection
    force_basis = R_a_inv(:, 1:n_af);
    [force_action, shape_margin] = forceControl(force_basis, ...
            eh_cones_goal_m, eh_cones_other_feasible_m);

    if ~isempty(force_action)
        disp('Force command:');
        disp(force_action);
        solution.eh_mode = eh_mode_goal;
        solution.n_af = n_af;
        solution.n_av = n_av;
        solution.w_av = w_av;
        solution.eta_af = -kForceMagnitude*force_action; % the minus sign comes from force balance
        solution.margin = min(shape_margin, margins(m));

        R_a_inv = R_a^-1;
        Cf_inv = vscale_inv*R_a_inv; Cf_inv = Cf_inv(:, 1:n_af);
        Cv_inv = vscale*R_a_inv; Cv_inv = Cv_inv(:, end-n_av+1:end);
        solution.R_a_inv = [Cf_inv Cv_inv];
        solution.R_a = solution.R_a_inv^-1;

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
        mode_text = printModes(solutions{i}.eh_mode, false);
        fprintf("Mode %d: %s\nmargin: %f\n", i, mode_text, solutions{i}.margin);
        margins(i) = solutions{i}.margin;
    end

    disp('Best solution:');
    [~, best_solution_id] = max(margins);
    solution = solutions{best_solution_id};
    V_T = solution.R_a_inv*[zeros(solution.n_af,1); solution.w_av];
    F_T = solution.R_a_inv*[solution.eta_af; zeros(solution.n_av, 1)];
    disp('R_a:');
    disp(solution.R_a);
    disp('V_T:');
    disp(V_T);
    disp('F_T:');
    disp(F_T);
end


return;
