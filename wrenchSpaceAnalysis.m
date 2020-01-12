function wrenchSpaceAnalysis(kFrictionE, kFrictionH, CP_e, CN_e, CP_O_h, CN_O_h, R_WO, p_WO, G, b_G, e_mode_goal, h_mode_goal)

kForceMagnitude = 5; % newton

Ne = size(CP_e, 2);
Nh = size(CP_O_h, 2);

R_OW = R_WO';
p_OW = -R_OW*p_WO;
Adj_OW = SE22Adj(R_OW(1:2,1:2), p_OW(1:2));

CP_h = CP_O_h;
CN_h = CN_O_h;
for i = 1:Nh
    CP_h(:, i) = R_WO*CP_O_h(:, i) + p_WO;
    CN_h(:, i) = R_WO*CN_O_h(:, i);
end

[Jac_e, Jac_h] = getWholeJacobian(CP_e, CN_e, CP_O_h, CN_O_h, Adj_OW);
[Jacf_e, Jacf_h] = getWholeJacobianFrictional(CP_e, CN_e, kFrictionE, CP_O_h, CN_O_h, kFrictionH, Adj_OW);

disp('Goal Mode:');
eh_mode_goal = [e_mode_goal; h_mode_goal];
goal_text = printModes(eh_mode_goal);
%%
%% Hybrid Servoing
%%

fprintf("###############################################\n");
fprintf("##             Hybrid Servoing               ##\n");
fprintf("###############################################\n");
%   velocity vector: v = [vo, vh]
% goal velocity
dims.Actualized = 3;
dims.UnActualized = 3;
dims.SlidingFriction = 0;

N_all = getJacobianFromContacts(e_mode_goal, h_mode_goal, Jac_e, Jac_h);

[n_av, n_af, R_a, w_av] = solvehfvc(dims, N_all, G, ...
    b_G, [], [], [], [], [], 3, false);

fprintf('The force controlled dimension is %d\n', n_af);
fprintf('The velocity controlled dimension is %d\n', n_av);

% make sure all velocity commands >= 0
R_id_flip = find(w_av < 0);
R_a(n_af + R_id_flip, :) = - R_a(n_af + R_id_flip, :);
w_av(R_id_flip) = - w_av(R_id_flip);

% % debug: show a crashing mode
% R_atemp = R_a;
% R_a(1, :) = R_atemp(3,:);
% R_a(3, :) = -R_atemp(1,:);
% R_a = R_a*aa2SO3(0.1, [0 0 1]);

Cv = [zeros(n_av, 3), R_a(n_af + 1: end, :)];
b_C = w_av;

R_a_inv = R_a^-1;

fprintf("###############################################\n");
fprintf("##              Mode Enumeration             ##\n");
fprintf("###############################################\n");

%
%% Use velocity command:
%  Filter out modes that are impossible:
%   1. If Nv = 0 doesn't have a solution, mark this mode as incompatible;
%       Otherwise, compute the solution.
%   2. If that solution satisfies Nu v >= 0, mark as feasible. Otherwise it is
%       impossible.

% Enumerate contact modes
e_modes = contact_mode_enumeration(CP_e(1:2,:), CN_e(1:2,:), true);
h_modes = contact_mode_enumeration(CP_O_h(1:2,:), CN_O_h(1:2,:), true);
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
disp(['Contact mode enumeration: ' num2str(size(e_modes, 2)*size(h_modes, 2)) ' modes.']);

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

% cone of possible actuation forces
%   force controlled direction can have both positive and negative forces
W_action = [R_a_inv(:, 1:n_af), -R_a_inv];

TOL = 1e-7;
feasible_mode_count = 0;
velocity_filtered_count = 0;
force_filtered_count = 0;
flag_crashing = false;

for i = 1:size(e_modes, 2)
    for j = 1:size(h_modes, 2)
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
                continue;
            end
            compatible = true;
        end

        [Je_, Jh_] = getFrictionalJacobianFromContacts(e_modes(:, i), h_modes(:, j), Jacf_e, Jacf_h);

        % check force balance
        R = coneIntersection(Je_', Jh_');
        if isempty(R) || norm(R) == 0
            force_filtered_count = force_filtered_count + 1;
            continue;
        end

        % check for feasible actuation
        R_ = coneIntersection(R, W_action);
        if isempty(R_) || norm(R_) == 0
            force_filtered_count = force_filtered_count + 1;
            continue;
        end

        % Crashing check:
        % Check if the infeasible mode contains v-controlled direction
        if ~compatible
            intersection = coneIntersection(R, W_action(:, end-n_av+1:end));
            if ~isempty(intersection) && norm(intersection) > TOL
                flag_crashing = true;
            end
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
        cone_generators{feasible_mode_count} = R_;
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
assert(n_af < 3);
assert(n_af > 0);
if flag_goal_is_feasible
    compatibilities(goal_id) = false;
    cone_generators_other_feasible = cone_generators(find(compatibilities));
    compatibilities(goal_id) = true;
    cone_generators_goal = cone_generators{goal_id};
    % project to force controlled plane/line
    force_basis = R_a_inv(:, 1:n_af);
    % pick a value
    projection_goal = force_basis'*cone_generators_goal;
    force_action = [];
    projection_goal_remains = projection_goal;

    force_action = mean(normc(projection_goal), 2);
    if ~isempty(cone_generators_other_feasible)
        if n_af == 1
            for i = 1:length(cone_generators_other_feasible)
                projection_other_feasible = force_basis'*cone_generators_other_feasible{i};
                if any(projection_other_feasible'*force_action > 0)
                    flag_force_region_feasible = false;
                    break;
                end
            end
        elseif n_af == 2
            for i = 1:length(cone_generators_other_feasible)
                projection_other_feasible = force_basis'*cone_generators_other_feasible{i};
                [~, projection_goal_remains] = coneIntersection2D( ...
                        projection_goal_remains, projection_other_feasible);
                if isempty(projection_goal_remains)
                    flag_force_region_feasible = false;
                    break;
                end
            end
            if flag_force_region_feasible
                force_action = mean(normc(projection_goal_remains), 2);
            end
        end
    end

    if flag_force_region_feasible
        disp('Force command:')
        disp(force_action);
    else
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
for i = 1:feasible_mode_count
    if i == goal_id
        fprintf('Mode %d: (Goal)\n', i);
    else
        fprintf('Mode %d:\n', i);
    end
    texts = printModes(eh_modes(:, i), false);
    texts(end) = [];
    if compatibilities(i) == true
        disp(texts);
        rcolor = [0.2*rand(), 0.5 + 0.5*rand(), 0.1+0.2*rand()];
        drawWrench(cone_generators{i}, rcolor, true, texts);

        cone_projection_i = force_basis*(force_basis')*cone_generators{i};
        drawWrench(cone_projection_i, rcolor, true);
    else
        disp([texts ' Infeasible']);
        drawWrench(cone_generators{i}, [0.2*rand(), 0.1+0.2*rand(), 0.5+0.5*rand()], true, texts);
    end
    fprintf('Margin: %f\n', margins(i));
    disp(cone_generators{i});
%     figure(i); clf(i); hold on;
%     drawWrench(cone_e{i},'r', true);
%     drawWrench(cone_h{i},'b', true);
end

if flag_goal_is_feasible
    % draw the goal region
    drawWrench(cone_generators{goal_id}, 'r', true);

    % draw the projections
    projection_goal_remains_3d = force_basis*projection_goal_remains;
    drawWrench(projection_goal_remains_3d, 'r', true);

    if flag_force_region_feasible
        faction = [[0;0;0] force_basis*force_action];
        plot3(faction(1,:), faction(2,:), faction(3,:), '-ro', 'markersize', 30);
    end
end


% todo: output action and score
n_av, n_af, R_a, w_av

