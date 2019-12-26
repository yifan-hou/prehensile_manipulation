function [Jac_e, Jac_h] = getWholeJacobian(CP_e, CN_e, CP_O_h, CN_O_h, p_WO, R_WO)

Ne = size(CP_e, 2);
Nh = size(CP_O_h, 2);

Jac_e = zeros(2*Ne, 6);
Jac_h = zeros(2*Nh, 6);

z = [0 0 1]';
for i = 1:Ne
    CN = CN_e(:,i)/norm(CN_e(:,i));
    pq = cross(CP_e(:,i), CN);
    CN_ = cross(z, CN);
    pq_ = cross(CP_e(:,i), CN_);
    Jac_e(2*i - 1, :) = [CN(1), CN(2), pq(3), 0, 0, 0];
    Jac_e(2*i,     :) = [CN_(1), CN_(2), pq_(3), 0, 0, 0];
end

adj_WO_inv_3 = SE22Adj(R_WO(1:2,1:2), p_WO(1:2));
for i = 1:Nh
    CN = CN_O_h(:,i)/norm(CN_O_h(:,i));
    pq = cross(CP_O_h(:,i), CN);
    c1 = [CN(1), CN(2), pq(3)]';
    c2 = c1'*adj_WO_inv_3;
    CN_ = cross(z, CN);
    pq_ = cross(CP_O_h(:,i), CN_);
    c1_ = [CN_(1), CN_(2), pq_(3)]';
    c2_ = c1_'*adj_WO_inv_3;
    Jac_h(2*i - 1, :) = [-c2, c2];
    Jac_h(2*i, :) = [-c2_, c2_];
end