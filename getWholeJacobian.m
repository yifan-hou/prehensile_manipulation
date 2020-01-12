function [Jac_e, Jac_h] = getWholeJacobian(CP_e, CN_e, CP_O_h, CN_O_h, adj_OW)

Ne = size(CP_e, 2);
Nh = size(CP_O_h, 2);

Jac_e = zeros(2*Ne, 6);
Jac_h = zeros(2*Nh, 6);

z = [0 0 1]';
for i = 1:Ne
    CN = CN_e(:,i)/norm(CN_e(:,i));
    pq = cross(CP_e(:,i), CN);
    CLeft = cross(z, CN);
    pq_ = cross(CP_e(:,i), CLeft);
    Jac_e(2*i - 1, :) = [CN(1), CN(2), pq(3), 0, 0, 0];
    Jac_e(2*i,     :) = [CLeft(1), CLeft(2), pq_(3), 0, 0, 0];
end

% CP_O gives constraints in object frame. But our hand velocity is in world frame.
% need transformation:
%   V_OH = V_OW + Adj_OW*V_WH = -Adj_OW*V_WO + Adj_OW*V_WH
% so
%   c*V_OH = (c*Adj_OW)*(V_WH - V_WO)
for i = 1:Nh
    CN = CN_O_h(:,i)/norm(CN_O_h(:,i));
    pq = cross(CP_O_h(:,i), CN);
    c1 = [CN(1), CN(2), pq(3)]';
    c2 = c1'*adj_OW;
    CLeft = cross(z, CN);
    pq_ = cross(CP_O_h(:,i), CLeft);
    c1_ = [CLeft(1), CLeft(2), pq_(3)]';
    c2_ = c1_'*adj_OW;
    Jac_h(2*i - 1, :) = [-c2, c2];
    Jac_h(2*i, :) = [-c2_, c2_];
end