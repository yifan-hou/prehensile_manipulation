function stabilityMarginOptimization(...
        kFrictionE, kFrictionH, CP_W_e, CN_W_e,...
        CP_H_h, CN_H_h, R_WH, p_WH, e_mode_goal, h_mode_goal)
addpath generated

tic
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

const_e = adj_WH';
const_h = -(adj_WH')*(adj_HW');

% X: a vector containing all the parameters.
%   X = [pe1x, pe1y, pe2x, pe2y, ...
%        ph1x, ph1y, ph2x, ph2y, ...]
Nx = (Ne + Nh)*2;
Nx_acc_Ne = Ne*2; % last position of Pe in X

NIter = 200;

margin_record = zeros(NIter, 1);
p_record = zeros(NIter, 1);
for iter = 1:NIter
    % =======================================
    %           compute the value
    % =======================================
    [Jacf_e, Jacf_h, NF_e, NF_h] = getWholeJacobianFrictional(CP_W_e, CN_W_e, kFrictionE, ...
            CP_H_h, CN_H_h, kFrictionH, adj_WH, adj_HW);

    % scale
    Jacf_e = Jacf_e * vscale;
    Jacf_h = Jacf_h * vscale;

    [Je_, Jh_, id_mode_e, id_mode_h] = getFrictionalJacobianFromContacts(...
        e_mode_goal, h_mode_goal, Jacf_e, Jacf_h);

    NJeRows = size(Je_, 1);
%     NJhRows = size(Jh_, 1);
    % check force balance
    R = coneIntersection(Je_', Jh_');
    if isempty(R) || norm(R) == 0
        disp('Error: force balance is lost.');
        return;
    end

    % compute cone stability margin
    % compute cone stability margin
    if sameConeCheck(Je_', R)
        [margin_, id_margin_R, id_margin_J1, id_margin_J2] = coneStabilityMargin(Jh_', R);
        id = 2;
    elseif sameConeCheck(Jh_', R)
        [margin_, id_margin_R, id_margin_J1, id_margin_J2] = coneStabilityMargin(Je_', R);
        id = 1;
    else
        [margin_e, id_margin_R_e, id_margin_Je1, id_margin_Je2] = coneStabilityMargin(Je_', R);
        [margin_h, id_margin_R_h, id_margin_Jh1, id_margin_Jh2] = coneStabilityMargin(Jh_', R);
        [margin_, id] = min([margin_e, margin_h]);
        if id == 1
            id_margin_R = id_margin_R_e;
            id_margin_J1 = id_margin_Je1;
            id_margin_J2 = id_margin_Je2;
        else
            id_margin_R = id_margin_R_h;
            id_margin_J1 = id_margin_Jh1;
            id_margin_J2 = id_margin_Jh2;
        end
    end


    fprintf('Margin: %f\n', margin_);
    margin_record(iter) = margin_;
    figure(1);clf(1);hold on;
    drawWrench(Je_','g', true);
    drawWrench(Jh_','b', true);
    view(-132, 26)
    drawnow

    % =======================================
    %           compute the gradient
    % =======================================
    %
    % step one: C - X
    %
    deriv_Ci_X = [];
    deriv_Cj_X = [];
    deriv_Ck_X = [];
    deriv_C1_X = [];
    deriv_C2_X = [];
    deriv_C3_X = [];
    deriv_C4_X = [];

    idm = idBackTrack_coneIntersection(Je_', Jh_', R);
    if id == 1
        % Ci, Cj belongs to cone_e
        Ci = normc(Je_(id_margin_J1, :)');
        Cj = normc(Je_(id_margin_J2, :)');
        deriv_Ci_X = getDerivCx(id_mode_e(id_margin_J1), NF_e, const_e, Nx, 0);
        deriv_Cj_X = getDerivCx(id_mode_e(id_margin_J2), NF_e, const_e, Nx, 0);
    else
        % Ci, Cj belongs to cone_h
        Ci = normc(Jh_(id_margin_J1, :)');
        Cj = normc(Jh_(id_margin_J2, :)');
        deriv_Ci_X = getDerivCx(id_mode_h(id_margin_J1), NF_h, const_h, Nx, Nx_acc_Ne);
        deriv_Cj_X = getDerivCx(id_mode_h(id_margin_J2), NF_h, const_h, Nx, Nx_acc_Ne);
    end
    % Ck
    CK_ide = idm(id_margin_R, 1:NJeRows);
    CK_idh = idm(id_margin_R, NJeRows+1:end);
    % Ck = [Je_(CK_ide, :);
    deriv_Cke_X = getDerivCx(id_mode_e(CK_ide), NF_e, const_e, Nx, 0);
    deriv_Ckh_X = getDerivCx(id_mode_h(CK_idh), NF_h, const_h, Nx, Nx_acc_Ne);
    Ck_is_vertex = sum(idm(id_margin_R, :)) == 1;
    if Ck_is_vertex
        % Ck is a vertex
        Ck = [Je_(CK_ide, :)' Jh_(CK_idh, :)']; % one of them is empty
        Ck = normc(Ck);
        deriv_Ck_X = [deriv_Cke_X deriv_Ckh_X]; % one of them is empty
    else
        % Ck is C1 C2 C3 C4
        C12 = normc(Je_(CK_ide, :)');
        C34 = normc(Jh_(CK_idh, :)');
        C1 = C12(:, 1);
        C2 = C12(:, 2);
        C3 = C34(:, 1);
        C4 = C34(:, 2);
        deriv_C1_X = deriv_Cke_X(:, :, 1);
        deriv_C2_X = deriv_Cke_X(:, :, 2);
        deriv_C3_X = deriv_Ckh_X(:, :, 1);
        deriv_C4_X = deriv_Ckh_X(:, :, 2);
    end

    %
    % step two: Phi - C
    %   four possibilities
    %
    sign_cross = 1;
    if Ck_is_vertex
        if isempty(deriv_Cj_X)
            % ik
            jac_phi_C = jac_phi_ik(Ci, Ck);
            C = [Ci Ck];
        else
            % ijk
            if cross(Ci, Cj)'*Ck < 0
                sign_cross = -1;
            end
            jac_phi_C = jac_phi_ijk(Ci, Cj, Ck);
            C = [Ci Cj Ck];
        end
    else
        C1234 = cross(cross(C1, C2), cross(C3, C4));
        if mean(C1234'*[C1 C2 C3 C4]) < 0
            C1234 = - C1234;
            sign_cross = -1;
        end
        if isempty(deriv_Cj_X)
            % i1234
            jac_phi_C = jac_phi_i1234(Ci, C1, C2, C3, C4);
            C = [Ci C1 C2 C3 C4];
        else
            % ij1234
            if cross(Ci, Cj)'*C1234 < 0
                sign_cross = -sign_cross;
            end
            jac_phi_C = jac_phi_ij1234(Ci, Cj, C1, C2, C3, C4);
            C = [Ci Cj C1 C2 C3 C4];
        end
    end
    jac_phi_C = sign_cross * jac_phi_C;

    % remove the components that scales C
    jac_phi_C = jac_phi_C - bsxfun(@times, dot(jac_phi_C, C), C);
    assert(norm(dot(jac_phi_C, C)) < 1e-7);

    % compute the final gradient
    if Ck_is_vertex
        if isempty(deriv_Cj_X)
            % ik
            jac_phi = jac_phi_C(:)'*[deriv_Ci_X;deriv_Ck_X];
        else
            % ijk
            jac_phi = jac_phi_C(:)'*[deriv_Ci_X;deriv_Cj_X;deriv_Ck_X];
        end
    else
        if isempty(deriv_Cj_X)
            % i1234
            jac_phi = jac_phi_C(:)'*[deriv_Ci_X;deriv_C1_X;deriv_C2_X;deriv_C3_X;deriv_C4_X];
        else
            % ij1234
            jac_phi = jac_phi_C(:)'*[deriv_Ci_X;deriv_Cj_X;deriv_C1_X;deriv_C2_X;deriv_C3_X;deriv_C4_X];
        end
    end
%     disp('Gradient: ');
%     disp(jac_phi);
    % =======================================
    %           Update parameter
    % =======================================
%     delta = 0.005;
%     CP_H_h(1, 1) = CP_H_h(1, 1) + delta*jac_phi(5);
%     CP_H_h(1, 2) = CP_H_h(1, 2) + delta*jac_phi(7);
%     disp(CP_H_h);

    delta = 0.002;
    CP_H_h(1) = CP_H_h(1) + delta*jac_phi(5);
    p_record(iter) = CP_H_h(1);
    disp(CP_H_h(1));
end
toc

figure(1);clf(1);hold on;
plot(margin_record,'.-b','markersize',2);
title('margin');
figure(2);clf(2);hold on;
plot(p_record,'.-b','markersize',2);

save generated/llf.mat margin_record p_record

end

% id_C: id of contact screw in whole Jacobian of Je or Jh
%   if id_C has one entry, output is 3xNx
%   if id_C has multiple entry, output is 3xNxxNC
function deriv = getDerivCx(id_C, NF, const, Nx, Nprevious)
    if isempty(id_C)
        deriv = [];
        return;
    end
    getPid = @(id) ceil(id/2);

    NC = length(id_C);
    if NC == 1
        NF = NF(:, id_C);
        deriv_Ci_p = const*[0 0;0 0; NF(2) -NF(1)];
        deriv = zeros(3, Nx);
        id_point = getPid(id_C); % id of contact point
        deriv(:, (Nprevious+2*id_point-1):(Nprevious+2*id_point)) = deriv_Ci_p;
    else
        NF = NF(:, id_C);
        deriv_Ci_p = zeros(3, 2, NC);
        deriv_Ci_p(3, 1, :) =  NF(2, :);
        deriv_Ci_p(3, 2, :) =  -NF(1, :);
        deriv_Ci_p = mat_ten(const, deriv_Ci_p);
        deriv = zeros(3, Nx, NC);
        id_point = getPid(id_C); % id of contact point
        for c = 1:NC
            deriv(:, (Nprevious+2*id_point(c)-1):(Nprevious+2*id_point(c)), c) = deriv_Ci_p(:, :, c);
        end
    end
end