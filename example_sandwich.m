%
%             -----------
%          h2|          | h1
%        --------------------
%        |        Y         |
%        |        ^         |
%     e2 |        | O       | e1
%    =============|---> X ===========
clc;clear;
%%
%% Problem definition
%%

% Parameters
kFrictionH = 1;
kFrictionE = 0.5;

% list of contact points and contact normals
p_e1 = [1; 0; 0];
p_e2 = [-1; 0; 0];
n_e1 = [0; 1; 0];
n_e2 = [0; 1; 0];

p_O_h1 = [0.5; 1; 0];
p_O_h2 = [-0.5; 1; 0];
n_O_h1 = [0; -1; 0];
n_O_h2 = [0; -1; 0];

R_WO = eye(3);
p_WO = zeros(3, 1);

p_h1 = R_WO*p_O_h1 + p_WO;
p_h2 = R_WO*p_O_h2 + p_WO;
n_h1 = R_WO*n_O_h1;
n_h2 = R_WO*n_O_h2;

CP_e = [p_e1, p_e2];
CN_e = [n_e1, n_e2];
CP_h = [p_h1, p_h2];
CN_h = [n_h1, n_h2];
CP_O_h = [p_O_h1, p_O_h2];
CN_O_h = [n_O_h1, n_O_h2];

Ne = 2;
Nh = 2;

[Jac_e, Jac_h] = getWholeJacobian(CP_e, CN_e, CP_O_h, CN_O_h, p_WO, R_WO);

%%
%% Pick a desired mode
%% 0, 1 or 2
%%

% user/planner specify the desired mode for hand
h_mode = [1; 1]; % ff
% planner determine the desired mode for environment
e_mode = [1; 0]; % fs

%%
%% Hybrid Servoing
%%

%   velocity vector: v = [vo, vh]
% goal velocity
dims.Actualized = 3;
dims.UnActualized = 3;
dims.SlidingFriction = 0;

N_all = getJacobianFromContacts(e_mode, h_mode, Jac_e, Jac_h);

G = [1 1 0 0 0 0];
b_G = 0.1;

[n_av, n_af, R_a, w_av] = solvehfvc(dims, N_all, G, ...
    b_G, [], [], [], [], [], 3);

C = [zeros(n_av, 3), R_a(n_af + 1: end, :)];
b_C = w_av;

%%
%% Use velocit command:
%%  Filter out modes
%%  compute sliding direction for each contact
%%
% Enumerate contact modes
e_modes = contact_mode_enumeration_nsd(CP_e(1:2,:), CN_e(1:2,:), true);
h_modes = contact_mode_enumeration_nsd(CP_O_h(1:2,:), CN_O_h(1:2,:), true);
disp(['Contact mode enumeration: ' num2str(size(e_modes, 2)*size(h_modes, 2)) ' modes.']);


eh_modes = [];
for i = 1:size(e_modes, 2)
    for j = 1:size(h_modes, 2)
        if all(e_modes(:, i)==0)
            continue;
        end
        if all(h_modes(:, j)==0)
            continue;
        end
        N = getJacobianFromContacts(e_modes(:, i), h_modes(:, j), Jac_e, Jac_h);
        % filter out modes using velocity command
        if rank([N;C]) - rank(N) > 0
            % compute possible sliding directions
            % these computations are based on 'Criteria for Maintaining Desired Contacts for Quasi-Static Systems'
            Lambda_bar = [C; N];
            b_Lambda_bar = [b_C; zeros(size(N,1), 1)];
            v_star = Lambda_bar\b_Lambda_bar;

            e_mode_ij = e_modes(:, i);
            h_mode_ij = h_modes(:, j);
            for kk = 1:length(e_mode_ij)
                if e_mode_ij(kk) == 2
                    Mx = Jac_e(2*kk, :);
                    v_Lx = Mx*v_star;
                    if v_Lx > 1e-10
                        e_mode_ij(kk) = 3;
                    elseif v_Lx < -1e-10
                        e_mode_ij(kk) = 2;
                    else
                        e_mode_ij(kk) = 1;
                    end
                end
            end
            for kk = 1:length(h_mode_ij)
                if h_mode_ij(kk) == 2
                    Mx = Jac_h(2*kk, :);
                    v_Lx = Mx*v_star;
                    if v_Lx > 1e-10
                        h_mode_ij(kk) = 3;
                    elseif v_Lx < -1e-10
                        h_mode_ij(kk) = 2;
                    else
                        h_mode_ij(kk) = 1;
                    end
                end
            end
            eh_modes = [eh_modes [e_mode_ij; h_mode_ij]];
        end
    end
end

eh_modes = unique(eh_modes','rows')';
disp(['Compatible with velocity command: ' num2str(size(eh_modes,2))]);

%% filter based on force balance
% friction cone approximation
[CPF_e, CNF_e] = frictionCone2D(CP_e, CN_e, kFrictionE);
[CPF_h, CNF_h] = frictionCone2D(CP_h, CN_h, kFrictionH);

% contact screws
W_h_all = contactScrew(CPF_h, CNF_h);
W_e_all = contactScrew(CPF_e, CNF_e);

eh_modes_feasible = [];

for i = 1:size(eh_modes, 2)
    W_e_i = [];
    W_h_i = [];
    for j = 1:Ne
        if eh_modes(j ,i) == 1 % fix
            W_e_i = [W_e_i W_e_all(:, 2*j-1:2*j)];
        elseif eh_modes(j ,i) == 2 % right
            W_e_i = [W_e_i W_e_all(:, 2*j)];
        elseif eh_modes(j, i) == 3 % left
            W_e_i = [W_e_i W_e_all(:, 2*j-1)];
        end
    end
    for j = 1:Nh
        if eh_modes(Ne+j ,i) == 1 % fix
            W_h_i = [W_h_i W_h_all(:, 2*j-1:2*j)];
        elseif eh_modes(Ne+j ,i) == 2 % right
            W_h_i = [W_h_i W_h_all(:, 2*j)];
        elseif eh_modes(Ne+j, i) == 3 % left
            W_h_i = [W_h_i W_h_all(:, 2*j-1)];
        end
    end

    % check force balance
    
end



%% Project the range of feasible modes onto the force controlled sub-space


%% find a force command, compute robustness score

%
% mode (sliding)
CPF_e_s = CPF_e(:, [1,3]);
CNF_e_s = CNF_e(:, [1,3]);
W_e_s = contactScrew(CPF_e_s, CNF_e_s);
drawWrench(-W_e_s([1 2 6], :), 'r');


% draw contact screws
figure(1); clf(1); hold on;
drawWrench(W_h([1 2 6], :), 'b');
drawWrench(-W_e([1 2 6], :), 'g');


% mode (pivoting)
CPF_e_s = CPF_e(:, [1,3]);
CNF_e_s = CNF_e(:, [1,3]);
W_e_s = contactScrew(CPF_e_s, CNF_e_s);
drawWrench(-W_e_s([1 2 6], :), 'r');


view(0, 90);
