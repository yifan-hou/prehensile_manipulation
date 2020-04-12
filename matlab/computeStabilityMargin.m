function margin = computeStabilityMargin(...
        kFrictionE, kFrictionH, CP_W_e, CN_W_e,...
        CP_H_h, CN_H_h, R_WH, p_WH, e_mode_goal, h_mode_goal)
% scaling for generalized velocity
% V = gvscale * V_scaled
kCharacteristicLength = 0.15;
vscale = diag([1 1 1/kCharacteristicLength]);

% N_ * V_scaled = 0, N*V = 0,
% -> N_ = N*gvscale

Ne = size(CP_W_e, 2);
Nh = size(CP_H_h, 2);

% add virtual z axis
assert(size(CP_W_e, 1) == 2);
CP_W_e = [CP_W_e; zeros(1, Ne)];
CN_W_e = [CN_W_e; zeros(1, Ne)];
CP_H_h = [CP_H_h; zeros(1, Nh)];
CN_H_h = [CN_H_h; zeros(1, Nh)];
R_WH = [R_WH [0; 0]; 0 0 1];
p_WH = [p_WH; 0];

R_HW = R_WH';
p_HW = -R_HW*p_WH;
adj_HW = SE22Adj(R_HW(1:2,1:2), p_HW(1:2));
adj_WH = SE22Adj(R_WH(1:2,1:2), p_WH(1:2));

[Jacf_e, Jacf_h] = getWholeJacobianFrictional(CP_W_e, CN_W_e, kFrictionE, ...
        CP_H_h, CN_H_h, kFrictionH, adj_WH, adj_HW);

% scale
Jacf_e = Jacf_e * vscale;
Jacf_h = Jacf_h * vscale;

[Je_, Jh_] = getFrictionalJacobianFromContacts(...
    e_mode_goal, h_mode_goal, Jacf_e, Jacf_h);

% check force balance
R = coneIntersection(Je_', Jh_');
if isempty(R) || norm(R) == 0
    margin = 0;
    return;
end

% compute cone stability margin
if sameConeCheck(Je_', R)
    margin = coneStabilityMargin(Jh_', R);
elseif sameConeCheck(Jh_', R)
    margin = coneStabilityMargin(Je_', R);
else
    margin_e = coneStabilityMargin(Je_', R);
    margin_h = coneStabilityMargin(Jh_', R);
    margin = min([margin_e, margin_h]);
end


% fprintf('Margin: %f\n', margin);