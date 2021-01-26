% J = [J_e, 0; -J_h, J_h]: 1 normal, kNumSlidingPlanes tangential; used by contact mode enumeration
% J_e = [N_e; T_e]
% J_h = [N_h; T_h]   N: normal, T: tangential (XY for 3D, left for 2D)
% Each contact contributes 1 normal, 2 (1 for planar problem) tangential constraints.
% eCone, hCone: each row is a wrench space generator created by an edge of a friction cone.
%   3D: Each contact contributes 2d + 1 edges; the last one is a copy of the first one
%   2D: Each contact contributes 2 edges (left, right)
% If planar, kNumSlidingPlanes must = 1, the computation puts everything on XY plane
%
% Left/right convention in 2d:
%    Object motion w.r.t. the environment
%
%    cone edge 1,      cone edge 2
%         -----\-------/------
%         |     \  N  /      |
%         |      \ ^ /       |
% T_e <---|       \|/        | ---> right sliding
%     =============|==================
%

function [N_e, T_e, N_h, T_h, eCone, eTCone, hCone, hTCone] = getWholeJacobian(CP_W_e, CN_W_e, ...
        CP_H_h, CN_H_h, adj_WH, adj_HW, kNumSlidingPlanes, kFrictionE, kFrictionH)

kDim = size(CP_W_e, 1);
assert((kDim == 2) || (kDim == 3));

if kDim == 3
    kWrenchDim = 6;
else
    kWrenchDim = 3;
    assert(kNumSlidingPlanes == 1);
end

Ne = size(CP_W_e, 2);
Nh = size(CP_H_h, 2);


N_e = zeros(Ne, kWrenchDim);
N_h = zeros(Nh, kWrenchDim);

if kDim == 3
    kEdgesPerContact = 2*kNumSlidingPlanes+1;
    T_e = zeros(Ne*2, kWrenchDim);
    T_h = zeros(Nh*2, kWrenchDim);
else
    kEdgesPerContact = 2;
    T_e = zeros(Ne, kWrenchDim);
    T_h = zeros(Nh, kWrenchDim);
end
eCone = zeros(kEdgesPerContact*Ne, kWrenchDim);
hCone = zeros(kEdgesPerContact*Nh, kWrenchDim);
eTCone = zeros(Ne*kNumSlidingPlanes, kWrenchDim);
hTCone = zeros(Nh*kNumSlidingPlanes, kWrenchDim);
z = [0 0 1]';

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
if kDim == 3
    SP = zeros(kNumSlidingPlanes, 3);
    for i = 1:kNumSlidingPlanes
        SP(i, 1:2) = [cos(pi/kNumSlidingPlanes*i), sin(pi/kNumSlidingPlanes*i)];
    end
    SP = [SP; -SP; SP(1, :)];
    % contact tangential directions in contact frame
    CT = zeros(kEdgesPerContact, 3);
    for i = 1:kEdgesPerContact
        CT(i, :) = cross(SP(i, :), z');
    end
end


if kDim == 3
    vr = rand(3, 1);
    vr = vr/norm(vr);
else
    vr = z;
end



mu_norm = sqrt(1 + kFrictionE^2);
for i = 1:Ne
    % contact normal
    CN = CN_W_e(:,i)/norm(CN_W_e(:,i));
    if kDim == 3
        c_Wei = [CN', cross(CP_W_e(:,i), CN)'];
    else
        c_Wei = [CN', cross2(CP_W_e(:,i)', CN')'];
    end
    cadj = c_Wei*adj_WH;
    N_e(i, :) = cadj;

    % contact tangential and friction cones
    if kDim == 3
        CX = cross(vr, CN); CX = CX/norm(CX);
        CY = cross(CN, CX); CY = CY/norm(CY);
        CT_e = CT(:, 1:2)*[CX'; CY']; % tangential directions
        CXY = [CX';CY'];
    else
        CX = cross(vr, [CN; 0]); CX = CX/norm(CX);
        CT_e = [CX'; -CX']; % left, right
        CT_e = CT_e(:, 1:2);
        CXY = CX(1:2)';
    end
    CCone = (kFrictionE*CT_e + ones(kEdgesPerContact, 1)*CN')/mu_norm; % friction cone edges

    CT_e = CT_e(1:kNumSlidingPlanes, :);
    if kDim == 3
        c_Wei = [CXY, cross([1;1]*CP_W_e(:,i)', CXY)];
        eTConei = [CT_e, cross(ones(kNumSlidingPlanes, 1)*CP_W_e(:,i)', CT_e)];
        eConei = [CCone, cross(ones(kEdgesPerContact, 1)*CP_W_e(:,i)', CCone)];
        T_e(2*(i-1)+1:2*i, :) = c_Wei*adj_WH;
    else
        c_Wei = [CXY, cross2(CP_W_e(:,i)', CXY)];
        eTConei = [CT_e, cross2(CP_W_e(:,i)', CT_e)];
        eConei = [CCone, cross2(ones(kEdgesPerContact, 1)*CP_W_e(:,i)', CCone)];
        T_e(i, :) = c_Wei*adj_WH;
    end

    eTCone(kNumSlidingPlanes*(i-1)+1:kNumSlidingPlanes*i, :) = eTConei*adj_WH;
    eCone((kEdgesPerContact)*(i-1)+1:(kEdgesPerContact)*i, :) = eConei*adj_WH;
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
    if kDim == 3
        c_Whi = [CN', cross(CP_H_h(:,i), CN)'];
    else
        c_Whi = [CN', cross2(CP_H_h(:,i)', CN')'];
    end
    cadj = c_Whi*adj_HW*adj_WH;
    N_h(i, :) = -cadj;

    % contact tangential and friction cones
    if kDim == 3
        CX = cross(vr, CN); CX = CX/norm(CX);
        CY = cross(CN, CX); CY = CY/norm(CY);
        CT_h = CT(:, 1:2)*[CX';CY']; % tangential directions
        CXY = [CX';CY'];
    else
        CX = cross(vr, [CN; 0]); CX = CX/norm(CX);
        CT_h = [CX'; -CX'];
        CT_h = CT_h(:, 1:2);
        CXY = CX(1:2)';
    end
    CCone = (kFrictionH*CT_h + ones(kEdgesPerContact, 1)*CN')/mu_norm; % friction cone edges

    CT_h = CT_h(1:kNumSlidingPlanes, :);
    if kDim == 3
        c_Whi = [CXY, cross([1;1]*CP_H_h(:,i)', CXY)];
        hTConei = [CT_h, cross(ones(kNumSlidingPlanes, 1)*CP_H_h(:,i)', CT_h)];
        hConei = [CCone, cross(ones(kEdgesPerContact, 1)*CP_H_h(:,i)', CCone)];
        T_h(2*(i-1)+1:2*i, :) = -c_Whi*adj_HW*adj_WH;
    else
        c_Whi = [CXY, cross2(CP_H_h(:,i)', CXY)];
        hTConei = [CT_h, cross2(CP_H_h(:,i)', CT_h)];
        hConei = [CCone, cross2(ones(kEdgesPerContact, 1)*CP_H_h(:,i)', CCone)];
        T_h(i, :) = -c_Whi*adj_HW*adj_WH;
    end
    
    hTCone(kNumSlidingPlanes*(i-1)+1:kNumSlidingPlanes*i, :) = -hTConei*adj_HW*adj_WH;
    hCone((kEdgesPerContact)*(i-1)+1:(kEdgesPerContact)*i, :) = -hConei*adj_HW*adj_WH;
end
