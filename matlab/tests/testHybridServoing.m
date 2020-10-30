clc;clear;
addpath ..
addpath ../algorithms


warning('off', 'MATLAB:rankDeficientMatrix');

% Parameters
NSamples = 100; %1000
num_seeds = 10;

% list of contact points and contact normals

% 0:separation 1:fixed 2/3: sliding
%          ne  nh  emodes  hmodes  ng  N
config_type = struct('ne', 1, 'nh', 1, 'ehmodes', [1 1], 'ng', 1);
configs2D = repmat(config_type, 4, 1);

configs2D(1).ne = 1; configs2D(1).nh = 1;
configs2D(1).ehmodes = [1 1;
                        2 1];
configs2D(1).ng = [2 2];

configs2D(2).ne = 1; configs2D(2).nh = 2;
configs2D(2).ehmodes = [1 1 1;
                        2 1 1];
configs2D(2).ng = [2 2];

configs2D(3).ne = 2; configs2D(3).nh = 1;
configs2D(3).ehmodes = [2 2 1];
configs2D(3).ng = [2];

configs2D(4).ne = 2; configs2D(4).nh = 2;
configs2D(4).ehmodes = [2 2 1 1];
configs2D(4).ng = [2];

kEMinX = 0;
kEMaxX = 0.2;
kEMinY = 0;
kEMaxY = 0.1;

kHMinX = 0;
kHMaxX = 0.2;
kHMinY = 0;
kHMaxY = 0.1;

kFrictionH = 0.7; % shouldn't matter
kFrictionE = 0.3; % shouldn't matter

Fg = [0 -10 0 0 0 0]';

Js = cell(20 * NSamples, 1);
Gs = cell(20 * NSamples, 1);
b_Gs = cell(20 * NSamples, 1);
As = cell(20 * NSamples, 1);
b_As = cell(20 * NSamples, 1);
count = 1;

for s = 1:size(configs2D, 1)
    disp(['Generating problems for setting ' num2str(s)]);
    ne = configs2D(s).ne;
    nh = configs2D(s).nh;
    for m = 1:size(configs2D(s).ehmodes, 1)
        emodes = configs2D(s).ehmodes(m, 1:ne)';
        hmodes = configs2D(s).ehmodes(m, ne+1:end)';
        ng = configs2D(s).ng(m);
        for n = 1:NSamples
            % randomly generate configs
            p_We = bsxfun(@plus, [kEMinX; kEMinY], diag([kEMaxX - kEMinX, kEMaxY - kEMinY]) * rand(2, ne));
            p_Hh = zeros(2, nh); % bsxfun(@plus, [kHMinX; kHMinY], diag([kHMaxX - kHMinX, kHMaxY - kHMinY]) * rand(2, nh));
            n_We = normalizeByCol([rand(1, ne) - 0.5; rand(1, ne)]);
            n_Hh = -normalizeByCol([rand(1, nh) - 0.5; rand(1, nh)]);

            angle = rand()*90-45; % deg
            R_WH = rotz_deg(angle);
            R_WH = R_WH(1:2, 1:2);
            p_WH = [0; kEMaxY];

            R_HW = R_WH';
            p_HW = -R_HW*p_WH;
            adj_HW = SE22Adj(R_HW, p_HW);
            adj_WH = SE22Adj(R_WH, p_WH);

            [N_e, T_e, N_h, T_h, eCone, eTCone, hCone, hTCone] = getWholeJacobian(p_We, n_We, ...
                p_Hh, n_Hh, adj_WH, adj_HW, 1, kFrictionE, kFrictionH);

            [N, Nu, normal_ids] = getJacobianFromContacts(emodes, hmodes, N_e, N_h, T_e, T_h);

            % goal
            G = rand(ng, 6);
            b_G = rand(ng, 1);
            % guard condition
            nLambda = size(N,1);
            A = eye(nLambda);
            A = -A(normal_ids, :);
            b_A = -5*ones(size(A,1),1);

            % save
            Js{count} = N;
            Gs{count} = G;
            b_Gs{count} = b_G;
            As{count} = A;
            b_As{count} = b_A;

            count = count + 1;
        end
    end
end

count = count - 1;


C1s = cell(count, 1);
C2s = cell(count, 1);

dims.Actualized = 3;
dims.UnActualized = 3;

time1.velocity = zeros(count,1);
time1.force = zeros(count,1);
time2.velocity = zeros(count, 1);
time2.force = zeros(count, 1);

number_of_optimal1 = 0;
number_of_optimal2 = 0;

COND_TOL = 0.001;
score1 = -ones(count,1);
score2 = -ones(count,1);

for p = 1:count
    disp([num2str(p) '/' num2str(count)]);
    J = Js{p};
    G = Gs{p};
    b_G = b_Gs{p};
    A = As{p};
    b_A = b_As{p};

    % ochs
    [action, time] = ochs(dims, J, G, b_G, Fg, A, b_A, J);
    C1s{p} = action.C;
    time1.velocity(p) = time.velocity*1000;
    time1.force(p) = time.force*1000;

    % HS
    [action, time] = hybrid_servoing(dims, J, G, b_G, Fg, A, b_A, num_seeds, J);
    C2s{p} = action.C;
%     M = blkdiag(zeros(3), eye(3));
%     U = null(J)';
%     U_bar = U*M;
%     U_bar = null(null(U_bar)')';
%     C2s{p} = U_bar;
    time2.velocity(p) = time.velocity*1000;
    time2.force(p) = time.force*1000;
    if ~isempty(C1s{p})
        % score1(p) = cond(null(J)'*normalizeByRow(C1s{p})');
        score1(p) = cond(normalizeByRow([J; C1s{p}]));

        if score1(p) - 1 < COND_TOL
            number_of_optimal1 = number_of_optimal1 + 1;
        end
    end

    if ~isempty(C2s{p})
        % score2(p) = cond(null(J)'*normalizeByRow(C2s{p})');
        score2(p) = cond(normalizeByRow([J; C2s{p}]));
        if score2(p) - 1 < COND_TOL
            number_of_optimal2 = number_of_optimal2 + 1;
        end
    end
end

disp('Score1  Score2');
disp([score1 score2]);
time1.velocity = time1.velocity(score1 > 0);
time1.force = time1.force(score1 > 0);
time2.velocity = time2.velocity(score2 > 0);
time2.force = time2.force(score2 > 0);

both_solved = (score1 > 0) & (score2 > 0);
% score1 = score1(score1 > 0);
% score2 = score2(score2 > 0);
number_of_solved1 = sum(score1 > 0);
number_of_solved2 = sum(score2 > 0);
average_cond1 = mean(score1(both_solved));
average_cond2 = mean(score2(both_solved));
best_cond1 = min(score1(both_solved));
best_cond2 = min(score2(both_solved));
worst_cond1 = max(score1(both_solved));
worst_cond2 = max(score2(both_solved));


average_vel_time1 = mean(time1.velocity);
average_vel_time2 = mean(time2.velocity);
best_vel_time1 = min(time1.velocity);
best_vel_time2 = min(time2.velocity);
worst_vel_time1 = max(time1.velocity);
worst_vel_time2 = max(time2.velocity);

average_force_time1 = mean(time1.force);
average_force_time2 = mean(time2.force);
best_force_time1 = min(time1.force);
best_force_time2 = min(time2.force);
worst_force_time1 = max(time1.force);
worst_force_time2 = max(time2.force);

average_time1 = average_force_time1 + average_vel_time1;
average_time2 = average_force_time2 + average_vel_time2;
best_time1 = min(time1.force + time1.velocity);
best_time2 = min(time2.force + time2.velocity);
worst_time1 = max(time1.force + time1.velocity);
worst_time2 = max(time2.force + time2.velocity);


disp('Total number of problems: ');
disp(count);
disp(['Solved problems 1: ' num2str(number_of_solved1) ', 2: ' num2str(number_of_solved2)]);
disp(['Optimal solutions 1: ' num2str(number_of_optimal1) ', 2: ' num2str(number_of_optimal2)]);
disp('Crashing Index:');
disp(['Average cond 1: ' num2str(average_cond1) ', 2: ' num2str(average_cond2)]);
disp(['best cond 1: ' num2str(best_cond1) ', 2: ' num2str(best_cond2)]);
disp(['worst cond 1: ' num2str(worst_cond1) ', 2: ' num2str(worst_cond2)]);
disp('Velocity Time');
disp(['Average velocity time 1: ' num2str(average_vel_time1) ', 2: ' num2str(average_vel_time2)]);
disp(['best velocity time 1: ' num2str(best_vel_time1) ', 2: ' num2str(best_vel_time2)]);
disp(['worst velocity time 1: ' num2str(worst_vel_time1) ', 2: ' num2str(worst_vel_time2)]);
disp('Force Time:')
disp(['Average force time 1: ' num2str(average_force_time1) ', 2: ' num2str(average_force_time2)]);
disp(['best force time 1: ' num2str(best_force_time1) ', 2: ' num2str(best_force_time2)]);
disp(['worst force time 1: ' num2str(worst_force_time1) ', 2: ' num2str(worst_force_time2)]);
disp('Speed Up')
disp(['Average Speedup velocity: ' num2str(average_vel_time2/average_vel_time1)]);
disp(['Average Speedup force: ' num2str(average_force_time2/average_force_time1)]);
disp(['Average Speedup overall: ' num2str((average_force_time2 + average_vel_time2)/(average_force_time1 + average_vel_time1))]);

