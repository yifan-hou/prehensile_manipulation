function solution = wrenchSpaceAnalysis(...
        kFrictionE, kFrictionH, CP_W_e, CN_W_e, ...
        CP_H_h, CN_H_h, R_WH, p_WH, G, b_G, e_mode_goal, h_mode_goal, kForceMagnitude)
solution = [];
TOL = 1e-7;

% scaling for generalized velocity
% V = gvscale * V_scaled
kCharacteristicLength = 0.10; % m
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

disp('Goal Mode:');
eh_mode_goal = [e_mode_goal; h_mode_goal];
printModes(eh_mode_goal);

%%
%% Hybrid Servoing
%%
fprintf("###############################################\n");
fprintf("##    Hybrid Servoing & Crashing Check       ##\n");
fprintf("###############################################\n");

[N_all, ~, Nue] = getJacobianFromContacts(e_mode_goal, h_mode_goal, Jac_e, Jac_h);
[n_av, n_af, R_a, R_a_inv, w_av,Cv, b_C] = hybridServoing(N_all, Nue, G, b_G);
if isempty(n_av)
    return;
end
W_action = [R_a_inv(:, 1:n_af), -R_a_inv];
% Crashing check
intersection = coneIntersection(R_all_f, W_action(:, end-n_av+1:end));
if ~isempty(intersection) && norm(intersection) > TOL
    disp('Crashing!');
    return;
end

fprintf("###############################################\n");
fprintf("##             Mode Filtering                ##\n");
fprintf("###############################################\n");

% How to filter out modes:
% 1. If NC degenerates, mark this mode as incompatible;
% 2. If nominal velocity under NC exists, and it cause the contact point to
%    slide in a different direction, remove this mode.
eh_modes = [];
compatibilities = [];
margins = [];

cone_generators = cell(size(e_modes, 2)*size(h_modes, 2), 1);
cone_e = cell(size(e_modes, 2)*size(h_modes, 2), 1);
cone_h = cell(size(e_modes, 2)*size(h_modes, 2), 1);


feasible_mode_count = 0;
velocity_filtered_count = 0;
force_filtered_count = 0;

for i = 1:size(e_modes, 2)
    for j = 1:size(h_modes, 2)
        printModes([e_modes(:, i); h_modes(:, j)]);
        % filter out modes using velocity command
        [N, Nu] = getJacobianFromContacts(e_modes(:, i), h_modes(:, j), Jac_e, Jac_h);
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
                velocity_filtered_count = velocity_filtered_count + 1;
                disp('Violates velocity inequality constraints. Discard');
                continue;
            end
            compatible = true;
        end

        [Je_, Jh_] = getFrictionalJacobianFromContacts(e_modes(:, i), h_modes(:, j), Jacf_e, Jacf_h);

        % check force balance
        R = coneIntersection(Je_', Jh_');
        if isempty(R) || norm(R) == 0
            force_filtered_count = force_filtered_count + 1;
            disp('Je, Jh has no intersections. Discard');
            continue;
        end

        % check cone stability margin
        margin_e = computeStabilityMargin(Je_', R);
        margin_h = computeStabilityMargin(Jh_', R);
        margin_ = min(margin_e, margin_h);

%         figure(1);clf(1);hold on;
%         printModes([e_modes(:, i); h_modes(:, j)]);
%         fprintf('Margin: %f\n', margin_);
%         drawWrench(Je_','g', true);
%         drawWrench(Jh_','b', true);
%         drawWrench(W_action,'k', true);

        feasible_mode_count = feasible_mode_count + 1;
        cone_generators{feasible_mode_count} = R;
        cone_e{feasible_mode_count} = Je_';
        cone_h{feasible_mode_count} = Jh_';
        eh_modes = [eh_modes [e_modes(:, i); h_modes(:, j)]];
        compatibilities = [compatibilities compatible];
        margins = [margins margin_];
    end
end

disp(['Velocity filtered: ' num2str(velocity_filtered_count)]);
disp(['Force filtered: ' num2str(force_filtered_count)]);
disp(['Remaining feasible modes: ' num2str(feasible_mode_count)]);

% trim
cone_generators = cone_generators(1:feasible_mode_count);
cone_e = cone_e(1:feasible_mode_count);
cone_h = cone_h(1:feasible_mode_count);

fprintf("###############################################\n");
fprintf("##                 Evaluation                ##\n");
fprintf("###############################################\n");
flag_goal_is_feasible = false;

goal_id = [];
for i = 1:feasible_mode_count
    mode_i = eh_modes(:, i);
    if all(mode_i == eh_mode_goal)
        goal_id = i;
        break;
    end
end
flag_goal_is_feasible = false;
if isempty(goal_id)
    disp('Goal mode is V-Impossible.');
elseif compatibilities(goal_id) == false
    disp('Goal mode is V-Infeasible.');
else
    flag_goal_is_feasible = true;
end

% action selection
flag_force_region_feasible = true;

if flag_goal_is_feasible
    force_basis = R_a_inv(:, 1:n_af);
    cone_generators_goal = cone_generators{goal_id};
    compatibilities(goal_id) = false;
    cone_generators_other_feasible = cone_generators(find(compatibilities));
    compatibilities(goal_id) = true;
    [force_action, shape_margin] = forceControl(force_basis, ...
            cone_generators_goal, cone_generators_other_feasible);
    
    if ~isempty(force_action)
        disp('Force command:')
        disp(force_action);
    else
        flag_force_region_feasible = false;
        disp('No distinguishable force command.');
    end
end

fprintf("###############################################\n");
fprintf("##                 Drawing                   ##\n");
fprintf("###############################################\n");

% exam all the modes
% draw contact screws
figure(1); clf(1); hold on;

% draw the force controlled subspace
drawWrench(W_action(:, 1:2*n_af), [0.8, 0.4, 0.4], false);
drawWrench(W_action(:, end-n_av+1:end), [0.4, 0.8, 0.8], true);

% draw the origin
plot3(0,0,0,'r*','markersize',35);

% draw cones
for i = 1:feasible_mode_count
    if i == goal_id
        fprintf('Mode %d: (Goal)\n', i);
        fprintf('Margin: %f\n', margins(i));
    else
        fprintf('Mode %d:\n', i);
    end
    texts = printModes(eh_modes(:, i), false);
    if compatibilities(i) == true
        disp(texts);
        rcolor = [0.2*rand(), 0.3 + 0.7*rand(), 0.1+0.4*rand()];
        drawWrench(cone_generators{i}, rcolor, true, texts);

        cone_projection_i = force_basis*(force_basis')*cone_generators{i};
        drawWrench(cone_projection_i, rcolor, true);
    else
        disp([texts ' Infeasible']);
%         drawWrench(cone_generators{i}, [0.2*rand(), 0.1+0.2*rand(), 0.5+0.5*rand()], true, texts);
    end

    disp(cone_generators{i});
%     figure(i); clf(i); hold on;
    % drawWrench(cone_e{i},'r', true);
    % drawWrench(cone_h{i},'b', true);
end


if flag_goal_is_feasible
    % draw the cone of goal mode
    drawWrench(cone_generators{goal_id}, 'r', true);

    if flag_force_region_feasible
        % draw the projections
        % projection_goal_remains_3d = force_basis*projection_goal_remains;
        % drawWrench(projection_goal_remains_3d, 'r', true);

        % draw the action
        faction = [[0;0;0] force_basis*force_action];
        plot3(faction(1,:), faction(2,:), faction(3,:), '-ro', 'markersize', 30);

        % record the solution
        solution.n_af = n_af;
        solution.n_av = n_av;
        solution.w_av = w_av;
        solution.eta_af = -kForceMagnitude*force_action; % the minus sign comes from force balance
        solution.margin = min(shape_margin, margins(goal_id));

        R_a_inv = R_a^-1;
        Cf_inv = vscale_inv*R_a_inv; Cf_inv = Cf_inv(:, 1:n_af);
        Cv_inv = vscale*R_a_inv; Cv_inv = Cv_inv(:, end-n_av+1:end);
        solution.R_a_inv = [Cf_inv Cv_inv];
        solution.R_a = solution.R_a_inv^-1;

        V_T = solution.R_a_inv*[zeros(n_af,1); solution.w_av];
        F_T = solution.R_a_inv*[solution.eta_af; zeros(n_av, 1)];
        disp('R_a:');
        disp(solution.R_a);
        disp('V_T:');
        disp(V_T);
        disp('F_T:');
        disp(F_T);
    end
end
