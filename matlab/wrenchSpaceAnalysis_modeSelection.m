
%!
%! Wrench space analysis.
%!
%!
%! @return     { description_of_the_return_value }
%!
function solution = wrenchSpaceAnalysis_modeSelection(...
        Normal_e, Normal_h, Tangent_e, Tangent_h, eCone_allFix, hCone_allFix, G, b_G, kForceMagnitude,...
        e_modes, h_modes, e_mode_goal, h_mode_goal)

solution = [];
TOL = 1e-7;

goal_mode_is_given = false;
if nargin == 15
    goal_mode_is_given = true;
end

% scaling for generalized velocity
% V = gvscale * V_scaled
kCharacteristicLength = 0.15;
vscale = diag([1 1 1/kCharacteristicLength]);
vscale_inv = diag([1 1 kCharacteristicLength]);
gvscale = diag([1 1 1/kCharacteristicLength 1 1 1/kCharacteristicLength]);
% N_ * V_scaled = 0, N*V = 0,
% -> N_ = N*gvscale
Normal_e = Normal_e * vscale;
Normal_h = Normal_h * vscale;
Tangent_e = Tangent_e * vscale;
Tangent_h = Tangent_h * vscale;
eCone_allFix = eCone_allFix * vscale;
hCone_allFix = hCone_allFix * vscale;
G = G*gvscale;

cone_allFix = coneIntersection(eCone_allFix', hCone_allFix');

fprintf("###############################################\n");
fprintf("##         Compute Stability Margins         ##\n");
fprintf("###############################################\n");
% EH cone intersection
%   * Compute safety margins
%   * get rid of infeasible cone
%

eh_modes = [];
margins = [];

eh_cones = cell(size(e_modes, 2)*size(h_modes, 2), 1);
Jacs = cell(size(e_modes, 2)*size(h_modes, 2), 1);
Jacus = cell(size(e_modes, 2)*size(h_modes, 2), 1);

eh_cone_feasible_mode_count = 0;
goal_id = 0;
for i = 1:size(e_modes, 2)
    for j = 1:size(h_modes, 2)
        [Je_, Jh_] = getFrictionalJacobianFromContacts(e_modes(:, i), h_modes(:, j), eCone_allFix, hCone_allFix);

        % check force balance
        R = coneIntersection(Je_', Jh_');
        if isempty(R) || norm(R) == 0
            continue;
        end

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
        if margin_ < TOL
            continue;
        end

        % compute velocity Jacobian
        [N, Nu] = getJacobianFromContacts(e_modes(:, i), h_modes(:, j), Normal_e, Normal_h, Tangent_e, Tangent_h);

        % figure(1);clf(1);hold on;
        % printModes([e_modes(:, i); h_modes(:, j)]);
        % fprintf('Margin: %f\n', margin_);
        % drawCone(Je_','g', true);
        % drawCone(Jh_','b', true);

        eh_modes = [eh_modes [e_modes(:, i); h_modes(:, j)]];
        margins = [margins margin_];

        eh_cone_feasible_mode_count = eh_cone_feasible_mode_count + 1;
        eh_cones{eh_cone_feasible_mode_count} = R;
        Jacs{eh_cone_feasible_mode_count} = N;
        Jacus{eh_cone_feasible_mode_count} = Nu;

        if goal_mode_is_given
            if (i == e_mode_goal) && (j == h_mode_goal)
                goal_id = eh_cone_feasible_mode_count;
            end
        end
    end
end

disp(['Modes with Margin > 0: ' num2str(eh_cone_feasible_mode_count)]);

if goal_mode_is_given && (goal_id == 0)
    disp("Failure: the Cone of the goal mode is empty.");
    return;
end

% trim
eh_cones = eh_cones(1:eh_cone_feasible_mode_count);
Jacs = Jacs(1:eh_cone_feasible_mode_count);
Jacus = Jacus(1:eh_cone_feasible_mode_count);

%%
%% Begin to check each cone
%%
fprintf("###############################################\n");
fprintf("##            Velocity Filtering             ##\n");
fprintf("###############################################\n");

solutions = cell(eh_cone_feasible_mode_count, 1);
solutions_count = 0;

for m = 1:eh_cone_feasible_mode_count
    if goal_mode_is_given
        m = goal_id;
    else
        disp('***********');
        fprintf("Mode %d of %d\n", m, eh_cone_feasible_mode_count);
        printModes(eh_modes(:, m));
    end

    eh_mode_goal = eh_modes(:, m);
    fprintf("====================================\n");
    fprintf("= Hybrid Servoing & Crashing Check =\n");
    fprintf("====================================\n");
    N_all = Jacs{m};
    [n_av, n_af, R_a, R_a_inv, w_av, Cv, b_C] = hybridServoing(N_all, G, b_G);
    if isempty(n_av)
        disp("Failure: Hybrid Servoing returns no solution.")
        if goal_mode_is_given
            return;
        else
            continue;
        end
    end

    V_control_directions = -R_a_inv(:, end-n_av+1:end);
    F_control_directions = [R_a_inv(:, 1:n_af) -R_a_inv(:, 1:n_af)];
    intersection = coneIntersection(cone_allFix, V_control_directions);
    % Crashing check
    if ~isempty(intersection) && norm(intersection) > TOL
        disp('Failure: Crashing. The mode is not feasible.');
        if goal_mode_is_given
            return;
        else
            continue;
        end
    end

    fprintf("=======================\n");
    fprintf("=== Mode Filtering ===\n");
    fprintf("=======================\n");
    % How to filter out modes:
    % 1. If NC degenerates, mark this mode as incompatible;
    % 2. If nominal velocity under NC exists, and it cause the contact point to
    %    slide in a different direction, remove this mode.
    feasibilities = false(eh_cone_feasible_mode_count, 1);

    % cone of possible actuation forces
    %   force controlled direction can have both positive and negative forces

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
        % drawCone(eh_cones{n},'g', true);
        % drawCone(W_action,'k', true);
        if compatible
            feasible_mode_count = feasible_mode_count + 1;
            feasibilities(n) = true;
        end
    end

    if flag_goal_infeasible
        disp('Goal mode violates velocity equality constraints. Discard this mode.');
        if goal_mode_is_given
            return;
        else
            continue;
        end
    end
    if flag_goal_impossible
        disp('Goal mode violates velocity inequality constraints. Discard this mode.');
        if goal_mode_is_given
            return;
        else
            continue;
        end
    end
    disp(['Remaining feasible modes: ' num2str(feasible_mode_count)]);

    drawWrenchSpace(cone_allFix, eCone_allFix, hCone_allFix, ...
            V_control_directions, F_control_directions, eh_cones, eh_modes, ...
            eh_cone_feasible_mode_count, feasibilities, goal_id);

    eh_cones_goal_m = eh_cones{m};
    feasibilities(m) = 0;
    eh_cones_other_feasible_m = eh_cones(feasibilities);

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

        if goal_mode_is_given
            break;
        end
    else
        disp('Failure: no distinguishable force command.');
        if goal_mode_is_given
            return;
        else
            continue;
        end
    end
end

fprintf("###############################################\n");
fprintf("##                  Results                  ##\n");
fprintf("###############################################\n");
fprintf("Total number of feasible modes found: %d\n", solutions_count);
solution = [];
if solutions_count > 0
    if goal_mode_is_given
        solution = solutions(1);
    else
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
    end

    printModes(solution.eh_mode);
    V_T = solution.R_a_inv*[zeros(solution.n_af,1); solution.w_av];
    F_T = solution.R_a_inv*[solution.eta_af; zeros(solution.n_av, 1)];
    disp('R_a:');
    disp(solution.R_a);
    disp('V_T:');
    disp(V_T);
    disp('F_T:');
    disp(F_T);
end
