% J = [J_e, 0; -J_h, J_h]: 1 normal, d tangential; used by contact mode enumeration
% J_e = [N_e; T_e]   N: normal, T: tangential
% J_h = [N_h; T_h]
% Each contact contributes 1 normal, d tangential constraints.
% eCone, hCone: each contact contributes 2d edges
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

eCone = zeros(2*kNumSlidingPlanes*Ne, 6);
hCone = zeros(2*kNumSlidingPlanes*Nh, 6);

% z = [0 0 1]';
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
    c_Wei = [CX', cross(CP_W_e(:,i), CX)';
             CY', cross(CP_W_e(:,i), CY)'];
    c_Wei = zeros(kNumSlidingPlanes, 6);
    eConei = zeros(2*kNumSlidingPlanes, 6);
    for j = 0:kNumSlidingPlanes-1
        Ct = cos(pi/kNumSlidingPlanes*j) * CX + sin(pi/kNumSlidingPlanes*j) * CY;
        c_Wei(j + 1, :) = [Ct', cross(CP_W_e(:,i), Ct)'];
        CconePlus = (kFrictionE*Ct + CN)/mu_norm;
        CconeNeg = (-kFrictionE*Ct + CN)/mu_norm;
        eConei(2*j + 1:2*j+2, :) = [CconePlus', cross(CP_W_e(:,i), CconePlus)';
                                    CconeNeg',  cross(CP_W_e(:,i), CconeNeg)'];
    end
    T_e(kNumSlidingPlanes*(i-1)+1:kNumSlidingPlanes*i, :) = c_Wei*adj_WH;
    eCone(2*kNumSlidingPlanes*(i-1)+1:2*kNumSlidingPlanes*i, :) = eConei*adj_WH;
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
    c_Wei = [CN', cross(CP_H_h(:,i), CN)'];
    cadj = c_Wei*adj_HW*adj_WH;
    N_h(i, :) = -cadj;

    % contact tangential
    CX = cross(vr, CN); CX = CX/norm(CX);
    CY = cross(CN, CX); CY = CY/norm(CY);
    c_Wei = [CX', cross(CP_H_h(:,i), CX)';
             CY', cross(CP_H_h(:,i), CY)'];
    cadj = c_Wei*adj_HW*adj_WH;
    c_Wei = zeros(kNumSlidingPlanes, 6);
    hConei = zeros(2*kNumSlidingPlanes, 6);
    for j = 0:kNumSlidingPlanes-1
        Ct = cos(pi/kNumSlidingPlanes*j) * CX + sin(pi/kNumSlidingPlanes*j) * CY;
        c_Wei(j + 1, :) = [Ct', cross(CP_H_h(:,i), Ct)'];
        CconePlus = (kFrictionH*Ct + CN)/mu_norm;
        CconeNeg = (-kFrictionH*Ct + CN)/mu_norm;
        hConei(2*j+1 : 2*j+2, :) = [CconePlus', cross(CP_H_h(:,i), CconePlus)';
                                    CconeNeg', cross(CP_H_h(:,i), CconeNeg)'];
    end
    T_h(kNumSlidingPlanes*(i-1)+1:kNumSlidingPlanes*i, :) = -c_Wei*adj_HW*adj_WH;
    hCone(2*kNumSlidingPlanes*(i-1)+1 : 2*kNumSlidingPlanes*i, :) = -hConei*adj_HW*adj_WH;
end

eCone = eCone';
hCone = hCone';