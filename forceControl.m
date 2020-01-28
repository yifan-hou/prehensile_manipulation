function [force_action, shape_margin] = forceControl(force_basis, ...
        eh_cones_goal, eh_cones_others)
projection_goal = force_basis'*eh_cones_goal;
force_action = mean(normc(projection_goal), 2);
shape_margin = inf;

if isempty(eh_cones_others)
    return;
end

projection_goal_remains = projection_goal;
projection_others = eh_cones_others;
for i = 1:length(eh_cones_others)
    projection_others{i} = force_basis'*eh_cones_others{i};
end

flag_force_region_feasible = true;
n_af = size(force_basis, 2);
if n_af == 1
    shape_margin = norm(force_action);
    for i = 1:length(eh_cones_others)
        if any(projection_others{i}'*force_action > 0)
            flag_force_region_feasible = false;
            break;
        end
    end
elseif n_af == 2
    % 1. subtract other regions
    for i = 1:length(eh_cones_others)
        [~, projection_goal_remains] = coneIntersection2D( ...
                projection_goal_remains, projection_others{i});
        if isempty(projection_goal_remains)
            flag_force_region_feasible = false;
            break;
        end
    end
    % 2. pick a force action in the remaining region
    if flag_force_region_feasible
        % pick the action that is far away from other regions.
        shape_margin = angBTVec(projection_goal_remains(:, 1), projection_goal_remains(:, 2));
        if length(eh_cones_others) > 1
            % just take the mean, stay away from both sides
            force_action = mean(normc(projection_goal_remains), 2);
            shape_margin = shape_margin/2;
        else
            % there is only one another cone.
            ang1 = angBTVec(projection_others{1}, projection_goal_remains(:, 1));
            ang2 = angBTVec(projection_others{1}, projection_goal_remains(:, 2));
            if ang1 > ang2
                force_action = normc(projection_goal_remains(:, 1));
            else
                force_action = normc(projection_goal_remains(:, 2));
            end
        end
    end
end

if flag_force_region_feasible == false
    force_action = [];
    shape_margin = [];
end
