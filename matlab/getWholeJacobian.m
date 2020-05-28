% J = [J_e, 0; -J_h, J_h]: 1 normal, d tangential; used by contact mode enumeration
% J_e = [N_e; T_e]   N: normal, T: tangential
% J_h = [N_h; T_h]
% Each contact contributes 1 normal, d tangential constraints.
% eCone, hCone: each row is a wrench space generator created by an edge of a friction cone.
% Each contact contributes 2d + 1 edges; the last one is a copy of the first one
% If kNumSlidingPlanes = 2, assume planar problem in XY plane.
function [N_e, T_e, N_h, T_h, eCone, hCone] = getWholeJacobian(CP_W_e, CN_W_e, ...
        CP_H_h, CN_H_h, adj_WH, adj_HW, kNumSlidingPlanes, kFrictionE, kFrictionH)

Ne = size(CP_W_e, 2);
Nh = size(CP_H_h, 2);

N_e = zeros(Ne, 6);
N_h = zeros(Nh, 6);
T_e = zeros(Ne*kNumSlidingPlanes, 6);
T_h = zeros(Nh*kNumSlidingPlanes, 6);

% kd = kNumSlidingPlanes + 1;
% J_e = zeros(kd*Ne, 6);
% J_h = zeros(kd*Nh, 6);

eCone = zeros((2*kNumSlidingPlanes+1)*Ne, 6);
hCone = zeros((2*kNumSlidingPlanes+1)*Nh, 6);

% The following code generates a list of sliding modes. This list controls the
% order of friction cone facets in eCone and hCone.
% For example, for kNumSlidingPlanes = 4, it should look like
% phase = [1 0 0 0
%          1 1 0 0
%          1 1 1 0
%          1 1 1 1
%          0 1 1 1
%          0 0 1 1
%          0 0 0 1
%          0 0 0 0];
% phase = zeros(2*kNumSlidingPlanes, kNumSlidingPlanes);
% for i = 1:2*kNumSlidingPlanes-1
%     id1_start = i - kNumSlidingPlanes + 1 ;
%     id1_end = i;
%     if id1_start < 1
%         id1_start = 1;
%     end
%     if id1_end > kNumSlidingPlanes
%         id1_end = kNumSlidingPlanes;
%     end
%     phase(i, id1_start:id1_end) = 1;
% end
%
% Compute the order from bit-wise code:
% if phase(1) == 1
%     id = sum(phase);
% else
%     id = 2*kNumSlidingPlanes - sum(phase);
% end

% sliding planes
SP = zeros(kNumSlidingPlanes, 3);
for i = 1:kNumSlidingPlanes
    SP(i, :) = [cos(pi/kNumSlidingPlanes*i), sin(pi/kNumSlidingPlanes*i), 0];
end
SP = [SP; -SP; SP(1, :)];
% contact tangential directions in contact frame
z = [0 0 1]';
CT = zeros(2*kNumSlidingPlanes+1, 3);
for i = 1:2*kNumSlidingPlanes+1
    CT(i, :) = cross(SP(i, :), z);
end

vr = rand(3, 1);
vr = vr/norm(vr);
mu_norm = sqrt(1 + kFrictionE^2);
for i = 1:Ne
    % contact normal
    CN = CN_W_e(:,i)/norm(CN_W_e(:,i));
    c_Wei = [CN', cross(CP_W_e(:,i), CN)'];
    cadj = c_Wei*adj_WH;
    N_e(i, :) = cadj;

    % contact tangential and friction cones
    CX = cross(vr, CN); CX = CX/norm(CX);
    CY = cross(CN, CX); CY = CY/norm(CY);
    CT_e = CT(:, 1:2)*[CX';CY']; % tangential directions
    CCone = (kFrictionE*CT_e + ones(2*kNumSlidingPlanes+1, 1)*CN')/mu_norm; % friction cone edges

    CT_e = CT_e(1:kNumSlidingPlanes, :);
    c_Wei = [CT_e, cross(ones(kNumSlidingPlanes, 1)*CP_W_e(:,i)', CT_e)];
    eConei = [CCone, cross(ones(2*kNumSlidingPlanes+1, 1)*CP_W_e(:,i)', CCone)];
    T_e(kNumSlidingPlanes*(i-1)+1:kNumSlidingPlanes*i, :) = c_Wei*adj_WH;
    eCone((2*kNumSlidingPlanes+1)*(i-1)+1:(2*kNumSlidingPlanes+1)*i, :) = eConei*adj_WH;
end

% CP_O gives constraints in object frame. But our hand velocity is in world frame.
% need transformation:
%   V_OH = V_OW + Adj_OW*V_WH = -Adj_OW*V_WO + Adj_OW*V_WH
% so
%   c*V_OH = (c*Adj_OW)*(V_WH - V_WO)
mu_norm = sqrt(1+kFrictionH^2);
for i = 1:Nh
    % contact normal
    CN = CN_H_h(:,i)/norm(CN_H_h(:,i));
    c_Whi = [CN', cross(CP_H_h(:,i), CN)'];
    cadj = c_Whi*adj_HW*adj_WH;
    N_h(i, :) = -cadj;

    % contact tangential and friction cones
    CX = cross(vr, CN); CX = CX/norm(CX);
    CY = cross(CN, CX); CY = CY/norm(CY);
    CT_h = CT(:, 1:2)*[CX';CY']; % tangential directions
    CCone = (kFrictionH*CT_h + ones(2*kNumSlidingPlanes+1, 1)*CN')/mu_norm; % friction cone edges

    CT_h = CT_h(1:kNumSlidingPlanes, :);
    c_Whi = [CT_h, cross(ones(kNumSlidingPlanes, 1)*CP_H_h(:,i)', CT_h)];
    hConei = [CCone, cross(ones(2*kNumSlidingPlanes+1, 1)*CP_H_h(:,i)', CCone)];
    T_h(kNumSlidingPlanes*(i-1)+1:kNumSlidingPlanes*i, :) = -c_Whi*adj_HW*adj_WH;
    hCone((2*kNumSlidingPlanes+1)*(i-1)+1:(2*kNumSlidingPlanes+1)*i, :) = -hConei*adj_HW*adj_WH;
end