clc;clear;
addpath ..
addpath ../algorithms

warning('off', 'MATLAB:rankDeficientMatrix');
warning('off', 'MATLAB:nearlySingularMatrix');

% Parameters
NSamples = 100;

num_seeds = [3, 10];

ill_threshold = 1000;

% list of contact points and contact normals

% 1:fixed 0: sliding
%          ne  nh  emodes  hmodes  ng  N
config_type = struct('ne', 1, 'nh', 1, 'ehmodes', [1 1], 'ng', 1);
configs3D = repmat(config_type, 7, 1);

configs3D(1).ne = 1; configs3D(1).nh = 1;
configs3D(1).ehmodes = [1 1;
                        0 1];
% configs3D(1).ng = [1 2];
configs3D(1).ng = [6 8];

configs3D(2).ne = 1; configs3D(2).nh = 2;
configs3D(2).ehmodes = [1 1 1;
                        0 1 1];
% configs3D(2).ng = [1 1];
configs3D(2).ng = [4 6];

configs3D(3).ne = 1; configs3D(3).nh = 3;
configs3D(3).ehmodes = [1 1 1 1;
                        0 1 1 1];
% configs3D(3).ng = [1 1];
configs3D(3).ng = [3 5];


configs3D(4).ne = 2; configs3D(4).nh = 1;
configs3D(4).ehmodes = [1 1 1;
                        1 0 1;
                        0 0 1];
% configs3D(4).ng = [1 2 2];
configs3D(4).ng = [4 5 7];

configs3D(5).ne = 2; configs3D(5).nh = 2;
configs3D(5).ehmodes = [1 1 1 1;
                        1 0 1 1;
                        0 0 1 1];
% configs3D(5).ng = [1 2 2];
configs3D(5).ng = [2 3 5];

configs3D(6).ne = 2; configs3D(6).nh = 3;
configs3D(6).ehmodes = [1 1 1 1 1;
                        1 0 1 1 1;
                        0 0 1 1 1];
% configs3D(6).ng = [1 2 4];
configs3D(6).ng = [1 2 4];

configs3D(7).ne = 3; configs3D(7).nh = 1;
configs3D(7).ehmodes = [1 1 0 1;
                        1 0 0 1;
                        0 0 0 1];
% configs3D(7).ng = [1 1 2];
configs3D(7).ng = [4 4 6];

configs3D(8).ne = 3; configs3D(8).nh = 2;
configs3D(8).ehmodes = [1 1 0 1 1;
                        1 0 0 1 1;
                        0 0 0 1 1];
% configs3D(8).ng = [1 1 2];
configs3D(8).ng = [2 2 4];

configs3D(9).ne = 3; configs3D(9).nh = 3;
configs3D(9).ehmodes = [1 1 0 1 1 1;
                        1 0 0 1 1 1;
                        0 0 0 1 1 1];
% configs3D(9).ng = [1 1 2];
configs3D(9).ng = [1 1 3];

kEMinX = 0;
kEMaxX = 0.2;
kEMinY = 0;
kEMaxY = 0.2;
kEMinZ = 0;
kEMaxZ = 0.1;

kHMinX = 0;
kHMaxX = 0.2;
kHMinY = 0;
kHMaxY = 0.2;
kHMinZ = 0;
kHMaxZ = 0.1;

obj_weight = 10;

kFrictionH = 0.7; % shouldn't matter
kFrictionE = 0.3; % shouldn't matter

Js = cell(90*NSamples,1);
Gs = cell(90*NSamples,1);
b_Gs = cell(90*NSamples,1);
As = cell(90 * NSamples,1);
b_As = cell(90 * NSamples,1);
count = 1;

% one finger
for s = 1:size(configs3D, 1)
    disp(['Generating problems for setting ' num2str(s)]);
    ne = configs3D(s).ne;
    nh = configs3D(s).nh;
    for m = 1:size(configs3D(s).ehmodes, 1)
        emodes = configs3D(s).ehmodes(m, 1:ne)';
        hmodes = configs3D(s).ehmodes(m, ne+1:end)';
        nRowsE = 3*sum(emodes == 1) + sum(emodes == 0);
        nRowsH = 3*sum(hmodes == 1) + sum(hmodes == 0);
        hmodes_combined = [hmodes; hmodes; hmodes];
        ng = configs3D(s).ng(m);
        for n = 1:NSamples
            % randomly generate configs
            p_We = bsxfun(@plus, [kEMinX; kEMinY; kEMinZ], diag([kEMaxX - kEMinX, kEMaxY - kEMinY, kEMaxZ - kEMinZ]) * rand(3, ne));
            n_We = normalizeByCol([rand(1, ne) - 0.5; rand(1, ne) - 0.5; rand(1, ne)]);

            p_Hh1 = bsxfun(@plus, [kHMinX; kHMinY; kHMinZ], diag([kHMaxX - kHMinX, kHMaxY - kHMinY, kHMaxZ - kHMinZ]) * rand(3, nh));
            n_Hh1 = -normalizeByCol([rand(1, nh) - 0.5; rand(1, nh) - 0.5; rand(1, nh)]);
            p_Hh2 = bsxfun(@plus, [kHMinX; kHMinY; kHMinZ], diag([kHMaxX - kHMinX, kHMaxY - kHMinY, kHMaxZ - kHMinZ]) * rand(3, nh));
            n_Hh2 = -normalizeByCol([rand(1, nh) - 0.5; rand(1, nh) - 0.5; rand(1, nh)]);
            p_Hh3 = bsxfun(@plus, [kHMinX; kHMinY; kHMinZ], diag([kHMaxX - kHMinX, kHMaxY - kHMinY, kHMaxZ - kHMinZ]) * rand(3, nh));
            n_Hh3 = -normalizeByCol([rand(1, nh) - 0.5; rand(1, nh) - 0.5; rand(1, nh)]);

            p_Hh = [p_Hh1 p_Hh2 p_Hh3];
            n_Hh = [n_Hh1 n_Hh2 n_Hh3];

            ax = rand(3,1);
            theta = rand()*45*pi/180;
            R_WH = aa2SO3(theta, ax);
            p_WH = [0; 0; kEMaxY];

            R_HW = R_WH';
            p_HW = -R_HW*p_WH;
            adj_HW = SE32Adj(R_HW, p_HW);
            adj_WH = SE32Adj(R_WH, p_WH);

            [N_e, T_e, N_h, T_h, eCone, eTCone, hCone, hTCone] = getWholeJacobian(p_We, n_We, ...
                p_Hh, n_Hh, adj_WH, adj_HW, 2, kFrictionE, kFrictionH);
            [J, normal_ids] = getJacobianFromContacts3D(emodes, hmodes_combined, N_e, N_h, T_e, T_h);


            J3 = zeros(nRowsE+3*nRowsH, 24);
            J3(:, 1:6) = J(:, 1:6);
            J3(nRowsE+1:nRowsE+nRowsH, 7:12) = J(nRowsE+1:nRowsE+nRowsH, 7:12);
            J3(nRowsE+nRowsH+1:nRowsE+2*nRowsH, 13:18) = J(nRowsE+nRowsH+1:nRowsE+2*nRowsH, 7:12);
            J3(nRowsE+2*nRowsH+1:end, 19:24) = J(nRowsE+2*nRowsH+1:end, 7:12);
            J1 = J3(1:nRowsE+nRowsH, 1:12);
            J2 = J3(1:nRowsE+2*nRowsH, 1:18);
            normal_ids1 = normal_ids(normal_ids < nRowsE + nRowsH);
            normal_ids2 = normal_ids(normal_ids < nRowsE + 2*nRowsH);
            normal_ids3 = normal_ids;

            % goal
            nullJ1 = null(J1)';
            nullJ2 = null(J2)';
            nullJ3 = null(J3)';
            G1 = (rand(ng, size(nullJ1, 1)) - 0.5) * nullJ1;
            G2 = (rand(ng, size(nullJ2, 1)) - 0.5) * nullJ2;
            G3 = (rand(ng, size(nullJ3, 1)) - 0.5) * nullJ3;
            b_G1 = rand(ng, 1) - 0.5;
            b_G2 = rand(ng, 1) - 0.5;
            b_G3 = rand(ng, 1) - 0.5;
            if rand() > 0.5
                G1(:,1:6) = 0;
                G2(:,1:6) = 0;
                G3(:,1:6) = 0;
            end
            % guard condition
            nLambda1 = size(J1,1);
            nLambda2 = size(J2,1);
            nLambda3 = size(J3,1);
            A1 = eye(nLambda1);
            A2 = eye(nLambda2);
            A3 = eye(nLambda3);
            A1 = -A1(normal_ids1, :);
            A2 = -A2(normal_ids2, :);
            A3 = -A3(normal_ids3, :);
            b_A1 = -5*ones(size(A1,1),1);
            b_A2 = -5*ones(size(A2,1),1);
            b_A3 = -5*ones(size(A3,1),1);

            % save
            Js{count} = J1;
            Gs{count} = G1;
            b_Gs{count} = b_G1;
            As{count} = A1;
            b_As{count} = b_A1;

%             Js{count} = J2;
%             Gs{count} = G2;
%             b_Gs{count} = b_G2;
%             As{count} = A2;
%             b_As{count} = b_A2;
%
%             Js{count} = J3;
%             Gs{count} = G3;
%             b_Gs{count} = b_G3;
%             As{count} = A3;
%             b_As{count} = b_A3;

%             count = count + 3;
            count = count + 1;
        end
    end
end

count = count - 1;

C1s = cell(count,1);
C2s = cell(count,1);
C3s = cell(count,1);
C4s = cell(count,1);

dims.Actualized = 6;
dims.UnActualized = 6;

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

score1 = -ones(count,1);
score2 = -ones(count,1);
score3 = -ones(count,1);
score4 = -ones(count,1);

for p = 1:count
%     disp([num2str(p) '/' num2str(count)]);
    J = Js{p};
    G = Gs{p};
    b_G = b_Gs{p};
    A = As{p};
    b_A = b_As{p};
    dims.Actualized = size(J,2) - 6;
    Fg = zeros(size(J,2),1);
    Fg(3) = -obj_weight;

    % solve
    [action, time] = ochs(dims, J, G, b_G, Fg, A, b_A, J);
    C1s{p} = action.C;
    time1.velocity(p) = time.velocity*1000;
    time1.force(p) = time.force*1000;

    [action, time] = ochs_max_v(dims, J, G, b_G, Fg, A, b_A, J);
    C2s{p} = action.C;
    time2.velocity(p) = time.velocity*1000;
    time2.force(p) = time.force*1000;

    [action, time] = hybrid_servoing(dims, J, G, b_G, Fg, A, b_A, num_seeds(1), J);
    C3s{p} = action.C;
    time3.velocity(p) = time.velocity*1000;
    time3.force(p) = time.force*1000;

    [action, time] = hybrid_servoing(dims, J, G, b_G, Fg, A, b_A, num_seeds(2), J);
    C4s{p} = action.C;
    time4.velocity(p) = time.velocity*1000;
    time4.force(p) = time.force*1000;

%     evaluation
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

