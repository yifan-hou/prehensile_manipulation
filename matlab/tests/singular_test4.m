% Point finger sliding test

clear;clc;
p_We = [0;
        0];
n_We = [0;
        1];

p_Hh = [0;
        0];

n_Hh = [0;
        -1];

R_WH = eye(2);
p_WH = [0; 0.3];

G = [1 0 0 0 0 0];
b_G = 0.1;

R_HW = R_WH';
p_HW = -R_HW*p_WH;
adj_HW = SE22Adj(R_HW, p_HW);
adj_WH = SE22Adj(R_WH, p_WH);

kFrictionE = 0.3;
kFrictionH = 0.8;
emodes = [1];
hmodes = [1];
[N_e, T_e, N_h, T_h, eCone, eTCone, hCone, hTCone] = getWholeJacobian(p_We, n_We, ...
    p_Hh, n_Hh, adj_WH, adj_HW, 1, kFrictionE, kFrictionH);
[J, ~, normal_ids] = getJacobianFromContacts(emodes, hmodes, N_e, N_h, T_e, T_h);

J = J(:,1:5);
Jreg = orth(J')';

kN = 50;
scores = zeros(kN,1);
for i = 0:kN-1
    theta = i/kN * pi - pi/2;
    x = cos(theta);
    y = sin(theta);
    C = [0 0 0 x y];

    score = cond(normalizeByRow([Jreg; C]));
    disp([num2str(theta*180/pi) ',    ' num2str(score)]);
    scores(i+1) = score;
end
disp('min: ');
disp(min(scores));
    