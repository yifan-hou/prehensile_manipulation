function [N_e, T_e, N_h, T_h, eCone_allFix, eTCone_allFix, hCone_allFix, hTCone_allFix, F_G] = preProcessing(...
        kFrictionE, kFrictionH, kNumSlidingPlanes, kObjWeight, CP_W_e, CN_W_e, CP_H_h, ...
        CN_H_h, R_WH, p_WH, CP_W_G)

R_HW = R_WH';
p_HW = -R_HW*p_WH;
if size(R_HW, 1) == 2
    adj_HW = SE22Adj(R_HW, p_HW);
    adj_WH = SE22Adj(R_WH, p_WH);
else
    adj_HW = SE32Adj(R_HW, p_HW);
    adj_WH = SE32Adj(R_WH, p_WH);
end

[N_e, T_e, N_h, T_h, eCone_allFix, eTCone_allFix, hCone_allFix, hTCone_allFix] = getWholeJacobian(CP_W_e, CN_W_e, ...
        CP_H_h, CN_H_h, adj_WH, adj_HW, kNumSlidingPlanes, kFrictionE, kFrictionH);

% gravity
z = [0 0 -1]';
W_G = [z; cross(CP_W_G, z)];
F_G = kObjWeight*(W_G')*adj_WH;
F_G = F_G'; % reshape to column vector