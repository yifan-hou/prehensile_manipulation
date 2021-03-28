clc;clear;
addpath ..
addpath ../algorithms


warning('off', 'MATLAB:rankDeficientMatrix');

% Parameters
NSamples = 1000; %1000
num_seeds = [3, 10];

dims.Actualized = 3;
dims.UnActualized = 3;

ill_threshold = 1000;

% list of contact points and contact normals

% 0:separation 1:fixed 2/3: sliding
%          ne  nh  emodes  hmodes  ng  N
config_type = struct('ne', 1, 'nh', 1, 'ehmodes', [1 1], 'ng', 1);
configs2D = repmat(config_type, 4, 1);

configs2D(1).ne = 1; configs2D(1).nh = 1;
configs2D(1).ehmodes = [1 1;
                        2 1];
% configs2D(1).ng = [2 2];
configs2D(1).ng = [2 3];

configs2D(2).ne = 1; configs2D(2).nh = 2;
configs2D(2).ehmodes = [1 1 1;
                        2 1 1];
% configs2D(2).ng = [1 2];
configs2D(2).ng = [1 2];

configs2D(3).ne = 2; configs2D(3).nh = 1;
configs2D(3).ehmodes = [2 2 1];
% configs2D(3).ng = [2];
configs2D(3).ng = [2];

configs2D(4).ne = 2; configs2D(4).nh = 2;
configs2D(4).ehmodes = [2 2 1 1];
% configs2D(4).ng = [2];
configs2D(4).ng = [1];

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
            p_Hh = bsxfun(@plus, [kHMinX; kHMinY], diag([kHMaxX - kHMinX, kHMaxY - kHMinY]) * rand(2, nh));
            p_Hh = bsxfun(@minus, p_Hh, mean(p_Hh,2));
            n_We = normalizeByCol([rand(1, ne) - 0.5; rand(1, ne)]);
            n_Hh = -normalizeByCol([rand(1, nh) - 0.5; rand(1, nh)]);

            angle = rand()*70-35; % deg
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

%             figure(1);clf(1); hold on;
%             v_We = [p_We p_We+0.1*n_We];
%             v_Hh = [p_Hh p_Hh+0.1*n_Hh];
%             v_Wh = R_WH*v_Hh + p_WH;
%             xH = [p_WH p_WH + R_WH(:,1)*0.05];
%             yH = [p_WH p_WH + R_WH(:,2)*0.05];
%             plot(v_We(1,:), v_We(2,:), '-r');
%             plot(p_We(1), p_We(2), 'r.', 'markersize',15);
%             plot(v_Wh(1,:), v_Wh(2,:), '-g');
%             plot(v_Wh(1,1), v_Wh(2,1), 'g.', 'markersize',15);
%             plot(xH(1,:), xH(2,:), '-*b');
%             plot(yH(1,:), yH(2,:), '-ob');
%             plot(0,0, '.r', 'markersize',15);
%             axis equal

%             % debug
%             [action, time] = ochs(dims, N, G, b_G, Fg, A, b_A, N);
%             C1 = action.C;
%             [action, time] = hybrid_servoing(dims, N, G, b_G, Fg, A, b_A, num_seeds, N);
%             C2 = action.C;
%             M = blkdiag(zeros(3), eye(3));
%             U = null(N)';
%             U_bar = U*M;
%             C3 = null(null(U_bar)')';
%
%             Nreg = rref(N);
%             Nreg = Nreg(1:rank(N),:);
%             if ~isempty(C1)
%                 score1 = cond(normalizeByRow([Nreg; C1]));
%                 score2 = cond(normalizeByRow([Nreg; C2]));
%                 score3 = cond(normalizeByRow([Nreg; C3]));
%             end

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
C3s = cell(count, 1);
C4s = cell(count, 1);

time1.velocity = zeros(count,1);
time1.force = zeros(count,1);
time2.velocity = zeros(count, 1);
time2.force = zeros(count, 1);
time3.velocity = zeros(count, 1);
time3.force = zeros(count, 1);
time4.velocity = zeros(count, 1);
time4.force = zeros(count, 1);

number_of_ill1 = 0;
number_of_ill2 = 0;
number_of_ill3 = 0;
number_of_ill4 = 0;

COND_TOL = 0.001;
score1 = -ones(count,1);
score2 = -ones(count,1);
score3 = -ones(count,1);
score4 = -ones(count,1);

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

    % ochs-max-v
    [action, time] = ochs_max_v(dims, J, G, b_G, Fg, A, b_A, J);
    C2s{p} = action.C;
    time2.velocity(p) = time.velocity*1000;
    time2.force(p) = time.force*1000;

    % HS
    [action, time] = hybrid_servoing(dims, J, G, b_G, Fg, A, b_A, num_seeds(1), J);
    C3s{p} = action.C;
    time3.velocity(p) = time.velocity*1000;
    time3.force(p) = time.force*1000;

    % HS
    [action, time] = hybrid_servoing(dims, J, G, b_G, Fg, A, b_A, num_seeds(2), J);
    C4s{p} = action.C;
    time4.velocity(p) = time.velocity*1000;
    time4.force(p) = time.force*1000;

    Jreg = orth(J')';

    if ~isempty(C1s{p})
        score1(p) = cond(normalizeByRow([Jreg; C1s{p}]));
        if score1(p) > ill_threshold
            number_of_ill1 = number_of_ill1 + 1;
        end
    end
    if ~isempty(C2s{p})
        score2(p) = cond(normalizeByRow([Jreg; C2s{p}]));
        if score2(p) > ill_threshold
            number_of_ill2 = number_of_ill2 + 1;
        end
    end
    if ~isempty(C3s{p})
        score3(p) = cond(normalizeByRow([Jreg; C3s{p}]));
        if score3(p) > ill_threshold
            number_of_ill3 = number_of_ill3 + 1;
        end
    end
    if ~isempty(C4s{p})
        score4(p) = cond(normalizeByRow([Jreg; C4s{p}]));
        if score4(p) > ill_threshold
            number_of_ill4 = number_of_ill4 + 1;
        elseif score4(p) < score1(p) - 0.1
            disp('nono');
            input('nono');
        end
    end
    disp(['Score1: ' num2str(score1(p)) ', Score2 ' num2str(score2(p)) ', Score3 ' num2str(score3(p)) ', Score4 ' num2str(score4(p))]);
end

time1.velocity = time1.velocity(score1 > 0);
time1.force = time1.force(score1 > 0);
time2.velocity = time2.velocity(score2 > 0);
time2.force = time2.force(score2 > 0);
time3.velocity = time3.velocity(score3 > 0);
time3.force = time3.force(score3 > 0);
time4.velocity = time4.velocity(score4 > 0);
time4.force = time4.force(score4 > 0);

solved1 = (score1 > 0) & (score1 < ill_threshold);
solved2 = (score2 > 0) & (score2 < ill_threshold);
solved3 = (score3 > 0) & (score3 < ill_threshold);
solved4 = (score4 > 0) & (score4 < ill_threshold);
ill_conditioned1 = score1 >= ill_threshold;
ill_conditioned2 = score2 >= ill_threshold;
ill_conditioned3 = score3 >= ill_threshold;
ill_conditioned4 = score4 >= ill_threshold;

all_solved = solved1 & solved2 & solved3 & solved4;

number_of_solved1 = sum(solved1);
number_of_solved2 = sum(solved2);
number_of_solved3 = sum(solved3);
number_of_solved4 = sum(solved4);

number_of_ill_conditioned1 = sum(ill_conditioned1);
number_of_ill_conditioned2 = sum(ill_conditioned2);
number_of_ill_conditioned3 = sum(ill_conditioned3);
number_of_ill_conditioned4 = sum(ill_conditioned4);

average_cond1 = mean(score1(all_solved));
average_cond2 = mean(score2(all_solved));
average_cond3 = mean(score3(all_solved));
average_cond4 = mean(score4(all_solved));

average_vel_time1 = mean(time1.velocity);
average_vel_time2 = mean(time2.velocity);
average_vel_time3 = mean(time3.velocity);
average_vel_time4 = mean(time4.velocity);
worst_vel_time1 = max(time1.velocity);
worst_vel_time2 = max(time2.velocity);
worst_vel_time3 = max(time3.velocity);
worst_vel_time4 = max(time4.velocity);

average_force_time1 = mean(time1.force);
average_force_time2 = mean(time2.force);
average_force_time3 = mean(time3.force);
average_force_time4 = mean(time4.force);
worst_force_time1 = max(time1.force);
worst_force_time2 = max(time2.force);
worst_force_time3 = max(time3.force);
worst_force_time4 = max(time4.force);

average_time1 = average_force_time1 + average_vel_time1;
average_time2 = average_force_time2 + average_vel_time2;
average_time3 = average_force_time3 + average_vel_time3;
average_time4 = average_force_time4 + average_vel_time4;
worst_time1 = max(time1.force + time1.velocity);
worst_time2 = max(time2.force + time2.velocity);
worst_time3 = max(time3.force + time3.velocity);
worst_time4 = max(time4.force + time4.velocity);

disp('Total number of problems: ');
disp(count);
disp(['Solved problems 1: ' num2str(number_of_solved1) ', 2: ' num2str(number_of_solved2) ', 3: ' num2str(number_of_solved3) ', 4: ' num2str(number_of_solved4)]);
disp('Crashing Index:');
disp(['Average cond 1: ' num2str(average_cond1) ', 2: ' num2str(average_cond2) ', 3: ' num2str(average_cond3) ', 4: ' num2str(average_cond4)]);
disp(['ill-conditioned: ' num2str(number_of_ill_conditioned1) ', ' num2str(number_of_ill_conditioned2) ', ' num2str(number_of_ill_conditioned3) ', ' num2str(number_of_ill_conditioned4)]);
disp('Velocity Time');
disp(['Average velocity time 1: ' num2str(average_vel_time1) ', 2: ' num2str(average_vel_time2) ', 3: ' num2str(average_vel_time3) ', 4: ' num2str(average_vel_time4)]);
disp(['worst velocity time 1: ' num2str(worst_vel_time1) ', 2: ' num2str(worst_vel_time2) ', 3: ' num2str(worst_vel_time3) ', 4: ' num2str(worst_vel_time4)]);
disp('Force Time:')
disp(['Average force time 1: ' num2str(average_force_time1) ', 2: ' num2str(average_force_time2) ', 3: ' num2str(average_force_time3) ', 4: ' num2str(average_force_time4)]);
disp(['worst force time 1: ' num2str(worst_force_time1) ', 2: ' num2str(worst_force_time2) ', 3: ' num2str(worst_force_time3) ', 4: ' num2str(worst_force_time4)]);
disp('Speed Up')
disp('Comparing with HS3')
disp(['Average Speedup velocity: ochs1: ' num2str(average_vel_time3/average_vel_time1) ', ochs2: ' num2str(average_vel_time3/average_vel_time2)]);
disp(['Average Speedup force: ochs1: ' num2str(average_force_time3/average_force_time1) ', ochs2: ' num2str(average_force_time3/average_force_time2)]);
disp(['Average Speedup overall: ochs1: ' num2str((average_force_time3 + average_vel_time3)/(average_force_time1 + average_vel_time1)) ', ochs2: ' num2str((average_force_time3 + average_vel_time3)/(average_force_time2 + average_vel_time2))]);
disp('Comparing with HS10')
disp(['Average Speedup velocity: ochs1: ' num2str(average_vel_time4/average_vel_time1) ', ochs2: ' num2str(average_vel_time4/average_vel_time2)]);
disp(['Average Speedup force: ochs1: ' num2str(average_force_time4/average_force_time1) ', ochs2: ' num2str(average_force_time4/average_force_time2)]);
disp(['Average Speedup overall: ochs1: ' num2str((average_force_time4 + average_vel_time4)/(average_force_time1 + average_vel_time1)) ', ochs2: ' num2str((average_force_time3 + average_vel_time3)/(average_force_time2 + average_vel_time2))]);

