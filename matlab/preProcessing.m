function [N_e, T_e, N_h, T_h, eCone_allFix, hCone_allFix, F_G] = preProcessing(...
        kFrictionE, kFrictionH, kNumSlidingPlanes, kObjWeight, CP_H_e, CN_H_e, CP_H_h, ...
        CN_H_h, CP_H_G, z_H)

[N_e, T_e, eCone_allFix] = getWholeJacobian_simple(CP_H_e, CN_H_e, kNumSlidingPlanes, kFrictionE);
[N_h, T_h, hCone_allFix] = getWholeJacobian_simple(CP_H_h, CN_H_h, kNumSlidingPlanes, kFrictionH);

% gravity
if (kNumSlidingPlanes == 1)
    H_G = [z_H; cross2(CP_H_G', z_H')];
else
    H_G = [z_H; cross(CP_H_G, z_H)];
end
F_G = kObjWeight*(H_G');
F_G = F_G'; % reshape to column vector
