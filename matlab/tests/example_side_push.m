clc;clear;
addpath ../
addpath ../plots/
addpath ../utilities/

% Parameters
kFrictionH = 0.7;
kFrictionE = 0.25;

%%
%% Geometrical Problem definition
%%

kW = 0.0435; % object width
kH = 0.0435; % object height

% center of gravity
p_W_G = [0; kH/2; 0];
% list of contact points and contact normals
p_W_e1 = [kW/2; 0; 0];
p_W_e2 = [-kW/2; 0; 0];
n_W_e1 = [0; 1; 0];
n_W_e2 = [0; 1; 0];

p_H_h1 = [-kW/2; 0; 0];
p_H_h2 = [-kW/2; -kH; 0];
n_H_h1 = [1; 0; 0];
n_H_h2 = [1; 0; 0];


CP_W_e = [p_W_e1, p_W_e2];
CN_W_e = [n_W_e1, n_W_e2];
CP_H_h = [p_H_h1, p_H_h2];
CN_H_h = [n_H_h1, n_H_h2];

angle = 0; % deg
R_WH = rotz(angle);
angle_rad = angle*pi/180;
p_WH = p_W_e2 + kW*[cos(angle_rad); sin(angle_rad); 0]/2 + kH*[-sin(angle_rad); cos(angle_rad); 0];

kNumSlidingPlanes = 4;
[J_e, J_h, eCone_allFix, hCone_allFix] = preProcessing(...
        kFrictionE, kFrictionH, kNumSlidingPlanes, CP_W_e, CN_W_e, CP_H_h, CN_H_h, R_WH, p_WH);

% mode enumeration
[e_modes, h_modes] = sharedGraspModeEnumeration(CP_W_e, CN_W_e, CP_H_h, CN_H_h);

e_mode = int8([1; 1]); % rr
h_mode = int8([1; 1]); % ff

kObjWeight = 10; % newton
kContactForce = 30;
margin = computeStabilityMargin(0.1+rand(), kFrictionH, ...
    CP_W_e(1:2,:), CN_W_e(1:2,:), CP_H_h(1:2,:), CN_H_h(1:2,:), p_W_G(1:2,:), ...
    R_WH(1:2,1:2), p_WH(1:2), e_mode, h_mode, kObjWeight, kContactForce);
disp('Minimal disturbance to break the mode:');
disp(margin);

