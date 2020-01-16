% N = [Jac_e; Jac_h]
function [Jac_e, Jac_h] = getWholeJacobian(CP_W_e, CN_W_e, ...
        CP_H_h, CN_H_h, adj_WH, adj_HW)

Ne = size(CP_W_e, 2);
Nh = size(CP_H_h, 2);

Jac_e = zeros(2*Ne, 6);
Jac_h = zeros(2*Nh, 6);

z = [0 0 1]';
for i = 1:Ne
    CN = CN_W_e(:,i)/norm(CN_W_e(:,i));
    pq = cross(CP_W_e(:,i), CN);
    c1 = [CN(1), CN(2), pq(3)];
    c2 = c1*adj_WH;
    CLeft = cross(z, CN);
    pq_ = cross(CP_W_e(:,i), CLeft);
    c1_ = [CLeft(1), CLeft(2), pq_(3)];
    c2_ = c1_*adj_WH;
    Jac_e(2*i - 1, :) = [c2, 0, 0, 0];
    Jac_e(2*i,     :) = [c2_, 0, 0, 0];
end

% CP_O gives constraints in object frame. But our hand velocity is in world frame.
% need transformation:
%   V_OH = V_OW + Adj_OW*V_WH = -Adj_OW*V_WO + Adj_OW*V_WH
% so
%   c*V_OH = (c*Adj_OW)*(V_WH - V_WO)
for i = 1:Nh
    CN = CN_H_h(:,i)/norm(CN_H_h(:,i));
    pq = cross(CP_H_h(:,i), CN);
    c1 = [CN(1), CN(2), pq(3)];
    c2 = c1*adj_HW;
    CLeft = cross(z, CN);
    pq_ = cross(CP_H_h(:,i), CLeft);
    c1_ = [CLeft(1), CLeft(2), pq_(3)]';
    c2_ = c1_'*adj_HW;
    Jac_h(2*i - 1, :) = [c2, -c2];
    Jac_h(2*i, :) = [c2_, -c2_];
end