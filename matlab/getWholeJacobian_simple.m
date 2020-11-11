% Nov.11 2020 I changed the notation of Jh. Now Je and Je can be computed using the same format
% J = [J_e, 0; J_h, -J_h]: 1 normal, kNumSlidingPlanes tangential; used by contact mode enumeration
% J_e = [N_e; T]
% J_h = [N_h; T_h]   N: normal, T: tangential (XY for 3D, left for 2D)
% Each contact contributes 1 normal, 2 (1 for planar problem) tangential constraints.
% Cone: each row is a wrench space generator created by an edge of a friction cone.
%   3D: Each contact contributes 2d + 1 edges; the last one is a copy of the first one
%   2D: Each contact contributes 2 edges (left, right)
% If planar, kNumSlidingPlanes must = 1, the computation puts everything on XY plane
%
% Left/right convention in 2d:
%    Object motion w.r.t. the environment
%
%    cone edge 1,      cone edge 2
%         -----\-------/------
%         |     \  N  /      |
%         |      \ ^ /       |
% T <---|       \|/        | ---> right sliding
%     =============|==================
%

function [N, T, Cone, TCone] = getWholeJacobian_simple(CP_H, CN_H, kNumSlidingPlanes, kFriction)

kDim = size(CP_H, 1);
assert((kDim == 2) || (kDim == 3));

if kDim == 3
    kWrenchDim = 6;
else
    kWrenchDim = 3;
    assert(kNumSlidingPlanes == 1);
end

kNumContacts = size(CP_H, 2);
N = zeros(kNumContacts, kWrenchDim);

if kDim == 3
    kEdgesPerContact = 2*kNumSlidingPlanes+1;
    T = zeros(kNumContacts*2, kWrenchDim);
else
    kEdgesPerContact = 2;
    T = zeros(kNumContacts, kWrenchDim);
end
Cone = zeros(kEdgesPerContact*kNumContacts, kWrenchDim);
TCone = zeros(kNumContacts*kNumSlidingPlanes, kWrenchDim);
z = [0 0 1]';

% The following code generates a list of sliding modes. This list controls the
% order of friction cone facets in eCone and hCone.
% For example, for kNumSlidingPlanes = 4, it should look like
% phase = [1 0 0 0
%          1 1 0 0
%          1 1 1 0
%          1 1 1 1
%          0 1 1 1
%          0 0 1 1
%          0 0 0 1
%          0 0 0 0];
% phase = zeros(2*kNumSlidingPlanes, kNumSlidingPlanes);
% for i = 1:2*kNumSlidingPlanes-1
%     id1_start = i - kNumSlidingPlanes + 1 ;
%     id1_end = i;
%     if id1_start < 1
%         id1_start = 1;
%     end
%     if id1_end > kNumSlidingPlanes
%         id1_end = kNumSlidingPlanes;
%     end
%     phase(i, id1_start:id1_end) = 1;
% end
%
% Compute the order from bit-wise code:
% if phase(1) == 1
%     id = sum(phase);
% else
%     id = 2*kNumSlidingPlanes - sum(phase);
% end

% sliding planes
if kDim == 3
    SP = zeros(kNumSlidingPlanes, 3);
    for i = 1:kNumSlidingPlanes
        SP(i, 1:2) = [cos(pi/kNumSlidingPlanes*i), sin(pi/kNumSlidingPlanes*i)];
    end
    SP = [SP; -SP; SP(1, :)];
    % contact tangential directions in contact frame
    CT = zeros(kEdgesPerContact, 3);
    for i = 1:kEdgesPerContact
        CT(i, :) = cross(SP(i, :), z');
    end
end

if kDim == 3
%     vr = rand(3, 1);
    vr = [0.1 -0.5 0.9]';
    vr = vr/norm(vr);
else
    vr = z;
end


mu_norm = sqrt(1 + kFriction^2);
for i = 1:kNumContacts
    % contact normal
    CN = CN_H(:,i)/norm(CN_H(:,i));
    if kDim == 3
        contact_screw_Hi = [CN', cross(CP_H(:,i), CN)'];
    else
        contact_screw_Hi = [CN', cross2(CP_H(:,i)', CN')'];
    end
    N(i, :) = contact_screw_Hi;

    % contact tangential and friction cones
    if kDim == 3
        CX = cross(vr, CN); CX = CX/norm(CX);
        CY = cross(CN, CX); CY = CY/norm(CY);
        CT_ = CT(:, 1:2)*[CX'; CY']; % tangential directions
        CXY = [CX';CY'];
    else
        CX = cross(vr, [CN; 0]); CX = CX/norm(CX);
        CT_ = [CX'; -CX']; % left, right
        CT_ = CT_(:, 1:2);
        CXY = CX(1:2)';
    end
    CCone = (kFriction*CT_ + ones(kEdgesPerContact, 1)*CN')/mu_norm; % friction cone edges

    CT_ = CT_(1:kNumSlidingPlanes, :);
    if kDim == 3
        contact_screw_Hi = [CXY, cross([1;1]*CP_H(:,i)', CXY)];
        TConei = [CT_, cross(ones(kNumSlidingPlanes, 1)*CP_H(:,i)', CT_)];
        Conei = [CCone, cross(ones(kEdgesPerContact, 1)*CP_H(:,i)', CCone)];
        T(2*(i-1)+1:2*i, :) = contact_screw_Hi;
    else
        contact_screw_Hi = [CXY, cross2(CP_H(:,i)', CXY)];
        TConei = [CT_, cross2(CP_H(:,i)', CT_)];
        Conei = [CCone, cross2(ones(kEdgesPerContact, 1)*CP_H(:,i)', CCone)];
        T(i, :) = contact_screw_Hi;
    end

    TCone(kNumSlidingPlanes*(i-1)+1:kNumSlidingPlanes*i, :) = TConei;
    Cone((kEdgesPerContact)*(i-1)+1:(kEdgesPerContact)*i, :) = Conei;
end