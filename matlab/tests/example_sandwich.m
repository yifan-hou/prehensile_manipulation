%             -----------
%          h2|          | h1
%        --------------------
%        |        Y         |
%        |        ^         |
%     e2 |        | O       | e1
%    =============|---> X ===========
clc;clear;
addpath ../
addpath ../plots/
addpath ../utilities/

% Parameters
kFrictionH = 0.7;
kFrictionE = 0.25;
kForceMagnitude = 15;

%%
%% Geometrical Problem definition
%%

kW = 0.0435; % object width
kH = 0.0435; % object height

% center of gravity
p_W_G = [0; kH/2];
% list of contact points and contact normals
p_W_e1 = [kW/2; 0];
p_W_e2 = [-kW/2; 0];
n_W_e1 = [0; 1];
n_W_e2 = [0; 1];

p_H_h1 = [kW/2; 0];
p_H_h2 = [-kW/2; 0];
n_H_h1 = [0; -1];
n_H_h2 = [0; -1];

p_H_h3 = [kW/2; -kH/4];
n_H_h3 = [0; -1];

CP_W_e = [p_W_e1, p_W_e2];
CN_W_e = [n_W_e1, n_W_e2];
CP_H_h = [p_H_h1, p_H_h2];
CN_H_h = [n_H_h1, n_H_h2];

angle = 0; % deg
R_WH = rotz(angle);
R_WH = R_WH(1:2, 1:2);
angle_rad = angle*pi/180;
p_WH = p_W_e2 + kW*[cos(angle_rad); sin(angle_rad)]/2 + kH*[-sin(angle_rad); cos(angle_rad)];

%%
%% Geometrical Pre-processing
%%
kNumSlidingPlanes = 1; % for 2D problem
[J_e, J_h, T_e, T_h, eCone_allFix, hCone_allFix] = preProcessing(...
        kFrictionE, kFrictionH, kNumSlidingPlanes, CP_W_e, CN_W_e, CP_H_h, CN_H_h, R_WH, p_WH);

% mode enumeration
[e_modes, h_modes] = sharedGraspModeEnumeration(CP_W_e, CN_W_e, CP_H_h, CN_H_h);

%%
%% Goal
%% 0, 1 or 2
%%

% Palm Pivot
G = [0 0 1 0 0 0];
b_G = [0.1];
e_mode = int8([0; 1]); % sf
h_mode = int8([1; 1]); % ff

% % Palm Slide
% G = [1 0 0 0 0 0];
% b_G = [0.1];
% e_mode = int8([3; 3]); % ss
% h_mode = int8([1; 1]); % ff

% % Finger Pivot
% G = [0 0 1 0 0 0];
% b_G = [0.1];
% e_mode = int8([0; 1]); % sf
% h_mode = int8([0; 1]); % sf

% % Palm Pivot slide
% adj_WH = SE32Adj(R_WH, p_WH);
% GO = [1 0 -p_W_e2(2)]*adj_WH;
% G = [GO 0 0 0;
%      0 0 1 0 0 0];
% b_G = [-0.1; 0.5];
% e_mode = int8([0; 3]); % sl
% h_mode = int8([1; 1]); % ff

% % Test three contacts
% G = [1 0 0 0 0 0];
% b_G = [0.1];
% e_mode = int8([3; 3]); % ss
% h_mode = int8([1; 1; 1]); % fff


% G = [adj_WH eye(3); 0 0 0 0 0 1];
% b_G = [0;0;0;-0.1];
% e_mode = int8(1); % f
% h_mode = int8(1); % f
% CP_W_e = CP_W_e(:, 2);
% CN_W_e = CN_W_e(:, 2);
% CP_H_h = CP_H_h(:, 1);
% CN_H_h = CN_H_h(:, 2);


tic
solution = wrenchSpaceAnalysis_modeSelection(J_e, J_h, T_e, T_h, eCone_allFix', ...
        hCone_allFix', G, b_G, kForceMagnitude, e_modes, h_modes, e_mode, h_mode);
toc

% p_H_h1 = [kW/4; 0];
% p_H_h2 = [-kW/4; 0];
% CP_H_h = [p_H_h1, p_H_h2];
%
% stabilityMarginOptimization(kFrictionE, kFrictionH, ...
%         CP_W_e(1:2,:), CN_W_e(1:2,:), CP_H_h(1:2,:), CN_H_h(1:2,:), R_WH(1:2,1:2), p_WH(1:2), e_mode, h_mode);
%
%
% kObjWeight = 10; % newton
% kContactForce = 30;
% tic
% for i = 1:10
%     margin = computeStabilityMargin(0.1+rand(), kFrictionH, ...
%         CP_W_e(1:2,:), CN_W_e(1:2,:), CP_H_h(1:2,:), CN_H_h(1:2,:), p_W_G(1:2,:), ...
%         R_WH(1:2,1:2), p_WH(1:2), e_mode, h_mode, kObjWeight, kContactForce);
% end
% toc
