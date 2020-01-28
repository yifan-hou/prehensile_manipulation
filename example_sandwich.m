%             -----------
%          h2|          | h1
%        --------------------
%        |        Y         |
%        |        ^         |
%     e2 |        | O       | e1
%    =============|---> X ===========
clc;clear;
%%
%% Problem definition
%%

% Parameters
kFrictionH = 1.5;
kFrictionE = 0.35;

kW = 0.0435; % object width
kH = 0.0435; % object height

% list of contact points and contact normals
p_W_e1 = [kW/2; 0];
p_W_e2 = [-kW/2; 0];
n_W_e1 = [0; 1];
n_W_e2 = [0; 1];

p_H_h1 = [kW/2; 0];
p_H_h2 = [-kW/2; 0];
n_H_h1 = [0; -1];
n_H_h2 = [0; -1];

CP_W_e = [p_W_e1, p_W_e2];
CN_W_e = [n_W_e1, n_W_e2];
CP_H_h = [p_H_h1, p_H_h2];
CN_H_h = [n_H_h1, n_H_h2];

angle = 0; % deg
R_WH = rotz(angle); R_WH = R_WH(1:2, 1:2);
angle_rad = angle*pi/180;
p_WH = p_W_e2 + kW*[cos(angle_rad); sin(angle_rad)]/2 + kH*[-sin(angle_rad); cos(angle_rad)];

%%
%% Goal
%% 0, 1 or 2
%%

% Palm Pivot
G = [0 0 1 0 0 0];
b_G = [0.1];
e_mode = int8([0; 1]); % sf
h_mode = int8([1; 1]); % ff

% % Finger Pivot
% G = [0 0 1 0 0 0];
% b_G = [0.1];
% e_mode = int8([0; 1]); % sf
% h_mode = int8([0; 1]); % sf

% Palm Pivot slide
G = [0 0 1 0 0 0;
     1 0 0 0 0 0];
b_G = [0.05; -1];
e_mode = int8([0; 3]); % sl
h_mode = int8([1; 1]); % ff

% adj_WH = SE22Adj(R_WH, p_WH);
% G = [adj_WH eye(3); 0 0 0 0 0 1];
% b_G = [0;0;0;-0.1];
% e_mode = int8(1); % f
% h_mode = int8(1); % f
% CP_W_e = CP_W_e(:, 2);
% CN_W_e = CN_W_e(:, 2);
% CP_H_h = CP_H_h(:, 1);
% CN_H_h = CN_H_h(:, 2);

solution = wrenchSpaceAnalysis(kFrictionE, kFrictionH, CP_W_e, CN_W_e, ...
        CP_H_h, CN_H_h, R_WH, p_WH, G, b_G, e_mode, h_mode, 15);

% solution = wrenchSpaceAnalysis_modeSelection(kFrictionE, kFrictionH, ...
%             CP_W_e, CN_W_e, CP_H_h, CN_H_h, R_WH, p_WH, G, b_G, 15);

% p_H_h1 = [kW/4; 0];
% p_H_h2 = [-kW/4; 0];
% CP_H_h = [p_H_h1, p_H_h2];
% 
% stabilityMarginOptimization(kFrictionE, kFrictionH, ...
%         CP_W_e, CN_W_e, CP_H_h, CN_H_h, R_WH, p_WH, e_mode, h_mode);