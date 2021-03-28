% clear;
% clc
% num_of_frames = 50; % config.num_of_frames;
% ang_max = 50 * pi/180; % config.rotation_angle_deg * pi/180;
% kFrictionE = 0.2;
% kFrictionH = 0.9;
% for fr = 0:num_of_frames
%     theta = ang_max*fr/num_of_frames;
%     st = sin(theta);
%     ct = cos(theta);
% 
%     % solve hfvc
%     dx = 0; %0.12;
%     dy = 0; %-0.24;
%     p_We = [-1 1; 0 0] + [dx dx; dy dy]; % each column is a variable
%     n_We = [0 0; 1 1];
%     p_Hh = [0; 0];
%     n_Hh = [-st; -ct];
% 
%     R_WH = rotz(theta*180/pi);
%     R_WH = R_WH(1:2,1:2);
%     p_WH = [0; 1] + [dx; dy];
% 
%     emodes = [2 2];
%     hmodes = [1];
%     
%     R_HW = R_WH';
%     p_HW = -R_HW*p_WH;
%     adj_HW = SE22Adj(R_HW, p_HW);
%     adj_WH = SE22Adj(R_WH, p_WH);
% 
%     [N_e, T_e, N_h, T_h, ~, ~, ~, ~] = getWholeJacobian(p_We, n_We, ...
%         p_Hh, n_Hh, adj_WH, adj_HW, 1, kFrictionE, kFrictionH);
% 
%     [N, ~] = getJacobianFromContacts(emodes, hmodes, N_e, N_h, T_e, T_h);
% 
%     C0 = [0 0 0 0 0 1;
%         0 0 0 ct -st 0];
%     
%     n_a = 3;
%     n_u = 3;
%     n   = 6;
%     
%     M = [zeros(n_u, n_a); eye(n_a)];
%     U = null(N)';
%     
%     U_bar = U*M;
%     U_hat = chol(U_bar'*U_bar + 1e-15*eye(n_a));
%     
%     TOL = 0.0001; % this tolerance has to be big.
%     id = [];
%     for i = 1:n_a
%         if (norm(U_hat(n_a-i + 1,:)) > TOL)
%             id = [id n_a-i+1];
%         end
%     end
%     U_hat = U_hat(id, :);
%     
%     n_av = size(U_hat,1);
%     n_af = n_a - n_av;
%     
%     % C_bar = (U_hat\eye(n_av))';
%     C_bar = lsqminnorm(U_hat, eye(n_av))';
%     C = [zeros(n_av, n_u) C_bar];
%     C2 = [zeros(n_av, n_u) U_bar];
%     C3 = [0 0 0 0 0 1;
%           0 0 0 1 0 0];
%     cond0 = cond(U*normalizeByRow(C0)');
%     cond1 = cond(U*normalizeByRow(C)');
%     cond2 = cond(U*normalizeByRow(C2)');
%     cond3 = cond(U*normalizeByRow(C3)');
% 
%     disp([cond0 cond1 cond2 cond3]);
% end


clc;
% % load test.mat
% 

p_We = [0; 0]; % each column is a variable
n_We = [1; 1];
p_Hh = [0; 0];
n_Hh = [0; -1];

R_WH = [1 0;0 1];
p_WH = [0; 1];

emodes = [1];
hmodes = [2];


R_HW = R_WH';
p_HW = -R_HW*p_WH;
adj_HW = SE22Adj(R_HW, p_HW);
adj_WH = SE22Adj(R_WH, p_WH);

[N_e, T_e, N_h, T_h, eCone, eTCone, hCone, hTCone] = getWholeJacobian(p_We, n_We, ...
    p_Hh, n_Hh, adj_WH, adj_HW, 1, 0.5, 0.9);

[N, Nu, normal_ids] = getJacobianFromContacts(emodes, hmodes, N_e, N_h, T_e, T_h);

% N = N(:,1:5);
n_a = 3;
n_u = 3;
n   = 6;

M = [zeros(n_u, n_a); eye(n_a)];
U = null(N)';

if isempty(U)
    % no null space. The environment already fully constrained the problem
    return;
end

U_bar = U*M;
U_hat = chol(U_bar'*U_bar + 1e-15*eye(n_a));

TOL = 0.0001; % this tolerance has to be big.
id = [];
for i = 1:n_a
    if (norm(U_hat(n_a-i + 1,:)) > TOL)
        id = [id n_a-i+1];
    end
end
U_hat = U_hat(id, :);


n_av = size(U_hat,1);
n_af = n_a - n_av;

% C_bar = (U_hat\eye(n_av))';
C_bar = lsqminnorm(U_hat, eye(n_av))';
C = [zeros(n_av, n_u) C_bar];
C2 = U*blkdiag(zeros(3), eye(3));
C2(normByRow(C2) < 1e-10, :) = [];

cond1 = cond(U*normalizeByRow(C)');
cond2 = cond(U*normalizeByRow(C2)');

disp([cond1 cond2]);