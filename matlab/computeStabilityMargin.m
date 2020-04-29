function margin = computeStabilityMargin(...
        kFrictionE, kFrictionH, CP_W_e, CN_W_e,...
        CP_H_h, CN_H_h, CP_W_G, R_WH, p_WH, e_mode_goal, h_mode_goal, kObjWeight, kContactForce)
% scaling for generalized velocity
% V = gvscale * V_scaled
kCharacteristicLength = 0.15;

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
CP_W_G = [CP_W_G; 0];
R_WH = [R_WH [0; 0]; 0 0 1];
p_WH = [p_WH; 0];

R_HW = R_WH';
p_HW = -R_HW*p_WH;
adj_HW = SE22Adj(R_HW(1:2,1:2), p_HW(1:2));
adj_WH = SE22Adj(R_WH(1:2,1:2), p_WH(1:2));

[Jacf_e, Jacf_h] = getWholeJacobianFrictional(CP_W_e, CN_W_e, kFrictionE, ...
        CP_H_h, CN_H_h, kFrictionH, adj_WH, adj_HW);

% gravity
W_G = contactScrew(CP_W_G, [0 -1 0]');
F_G = kObjWeight*v3t2(W_G)'*adj_WH;


% scale torque to force
vscale = diag([1 1 1/kCharacteristicLength]);
Jacf_e = Jacf_e * vscale;
Jacf_h = Jacf_h * vscale;

[Je_, Jh_] = getFrictionalJacobianFromContacts(...
    e_mode_goal, h_mode_goal, Jacf_e, Jacf_h);

% scale kContactForce
Je_ = kContactForce * Je_;
Jh_ = - kContactForce * Jh_;

% compute minkowski sum
polytope = [0 0 0; Je_(1, :)];
for i = 2:size(Je_, 1)
    polytope = [polytope; bsxfun(@plus, polytope, Je_(i, :))];
end
% id = unique(convhulln(polytope), 'stable');
% polytope = polytope(id, :);
polytope = bsxfun(@plus, polytope, F_G);
for i = 1:size(Jh_, 1)
    polytope = [polytope; bsxfun(@plus, polytope, Jh_(i, :))];
end
id = unique(convhulln(polytope), 'stable');
polytope = polytope(id, :);

% polytope_e = Polyhedron('V', zeros(1, 3));
% for i = 1:size(Je_, 1)
%     V = [0 0 0; Je_(i, :)] + [1;1] * F_G;
%     polytope_e = polytope_e + Polyhedron('V', V);
% end
% polytope_h = Polyhedron('V', zeros(1, 3));
% for i = 1:size(Jh_, 1)
%     polytope_h = polytope_h + Polyhedron('V', [0 0 0; -Jh_(i, :)]);
% end
% polytope = polytope_e + polytope_h;

% find the shortest distance from the origin to the face of the polytope
polytope = Polyhedron('V', polytope);
polytope.computeHRep;
A = polytope.A;
b = polytope.b;
if any(b <= 1e-5)
    margin = 0;
    return;
end

margin = min(b./normByRow(A));

