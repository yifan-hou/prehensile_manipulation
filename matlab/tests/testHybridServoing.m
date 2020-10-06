clc;clear;
addpath ..
addpath ../algorithms

% Parameters
NSamples = 100;
num_seeds = 5;

% list of contact points and contact normals

% 0:separation 1:fixed 2/3: sliding
%          ne  nh  emodes  hmodes  ng  N
config_type = struct('ne', 1, 'nh', 1, 'ehmodes', [1 1], 'ng', 1);
configs2D = repmat(config_type, 4, 1);

configs2D(1).ne = 1; configs2D(1).nh = 1;
configs2D(1).ehmodes = [1 1;
                        2 1;
                        3 1];
configs2D(1).ng = [1 2 1];

configs2D(2).ne = 1; configs2D(2).nh = 2;
configs2D(2).ehmodes = [1 1 1;
                        2 1 1];
configs2D(2).ng = [1 2];

configs2D(3).ne = 2; configs2D(3).nh = 1;
configs2D(3).ehmodes = [2 2 1;
                        3 3 1];
configs2D(3).ng = [1 1];

configs2D(4).ne = 2; configs2D(4).nh = 2;
configs2D(4).ehmodes = [2 2 1 1;
                        3 3 1 1];
configs2D(4).ng = [1 1];

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


Js = cell(20 * NSamples);
Gs = cell(20 * NSamples);
b_Gs = cell(20 * NSamples);
count = 1;

for s = 1:size(configs2D, 1)
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
            n_We = normalizeByCol([rand(1, ne) - 0.5; rand(1, ne)]);
            n_Hh = -normalizeByCol([rand(1, nh) - 0.5; rand(1, nh)]);

            angle = rand()*90-45; % deg
            R_WH = rotz(angle);
            R_WH = R_WH(1:2, 1:2);
            p_WH = [0; kEMaxY];

            R_HW = R_WH';
            p_HW = -R_HW*p_WH;
            adj_HW = SE22Adj(R_HW, p_HW);
            adj_WH = SE22Adj(R_WH, p_WH);

            [N_e, T_e, N_h, T_h, eCone, eTCone, hCone, hTCone] = getWholeJacobian(p_We, n_We, ...
                p_Hh, n_Hh, adj_WH, adj_HW, 1, kFrictionE, kFrictionH);

            [N, Nu] = getJacobianFromContacts(emodes, hmodes, N_e, N_h, T_e, T_h);

            % goal
            G = rand(ng, 6);
            b_G = rand(ng, 1);
            % save

            Js{count} = N;
            Gs{count} = G;
            b_Gs{count} = b_G;
            count = count + 1;
        end
    end
end

count = count - 1;
disp('Total number of problems: ');
disp(count);

C1s = cell(count);
C2s = cell(count);

dims.Actualized = 3;
dims.UnActualized = 3;

% 1. OCHS
tic;
for p = 1:count
    J = Js{p};
    G = Gs{p};
    b_G = b_Gs{p};

    % ochs
    [C, b_C] = ochs(dims, J, G, b_G);
    C1s{p} = C;
end
time1 = toc;

% 2. old hybrid servoing, seeds = 3
tic;
for p = 1:count
    J = Js{p};
    G = Gs{p};
    b_G = b_Gs{p};

    [C, b_C] = hybrid_servoing(dims, J, G, b_G, num_seeds);
    C2s{p} = C;
end
time2 = toc;

%% Evaluation
number_of_solved1 = 0;
number_of_optimal1 = 0;
average_cond_of_solved_problems1 = 0;
number_of_solved2 = 0;
number_of_optimal2 = 0;
average_cond_of_solved_problems2 = 0;

COND_TOL = 0.001;
for p = 1:count
    J = Js{p};
    G = Gs{p};
    b_G = b_Gs{p};

    score1 = -1;
    if ~isempty(C1s{p})
        score1 = cond(null(J)'*C1s{p}');
        number_of_solved1 = number_of_solved1 + 1;
        if score1 - 1 < COND_TOL
            number_of_optimal1 = number_of_optimal1 + 1;
        end
        average_cond_of_solved_problems1 = average_cond_of_solved_problems1 + score1;
    end


    score2 = -1;
    if ~isempty(C2s{p})
        score2 = cond(null(J)'*C2s{p}');
        number_of_solved2 = number_of_solved2 + 1;
        if score2 - 1 < COND_TOL
            number_of_optimal2 = number_of_optimal2 + 1;
        end
        average_cond_of_solved_problems2 = average_cond_of_solved_problems2 + score2;
    end

    disp(['score 1: ' num2str(score1) ', score2: ' num2str(score2)]);
end

disp(['time 1: ' num2str(time1*1000) ', time2: ' num2str(time2*1000)]);
disp(['Speedup: ' num2str(time2/time1)]);
disp(['Solved problems 1: ' num2str(number_of_solved1) ', 2: ' num2str(number_of_solved2)]);
disp(['Optimal solutions 1: ' num2str(number_of_optimal1) ', 2: ' num2str(number_of_optimal2)]);
disp(['Aver cond 1: ' num2str(average_cond_of_solved_problems1/count) ', 2: ' num2str(average_cond_of_solved_problems2/count)]);
