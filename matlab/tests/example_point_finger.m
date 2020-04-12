%
%                 --
%                 h1
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
kFrictionH = 1;
kFrictionE = 0.3;

kW = 0.1; % object width
kH = 0.1; % object height

% list of contact points and contact normals
p_W_e1 = [kW/2; 0];
p_W_e2 = [-kW/2; 0];
n_W_e1 = [0; 1];
n_W_e2 = [0; 1];

p_H_h1 = [0; 0];
n_H_h1 = [0; -1];

CP_W_e = [p_W_e1, p_W_e2];
CN_W_e = [n_W_e1, n_W_e2];
CP_H_h = [p_H_h1];
CN_H_h = [n_H_h1];


R_WH = eye(2);
p_WH = [0; kH];

%%
%% Goal
%% 0, 1 or 2

% % sff
% G = [0 0 1 0 0 0];
% b_G = 0.1;
% e_mode = int8([0; 1]); % fs
% h_mode = int8(1); % f

% llf
G = [1 0 0 0 0 0];
b_G = -0.1;
e_mode = int8([3; 3]); % ll
h_mode = int8(1); % f


% wrenchSpaceAnalysis(kFrictionE, kFrictionH, CP_W_e, CN_W_e, ...
%         CP_H_h, CN_H_h, R_WH, p_WH, G, b_G, e_mode, h_mode);

solution = wrenchSpaceAnalysis_modeSelection(kFrictionE, kFrictionH, ...
        CP_W_e, CN_W_e, CP_H_h, CN_H_h, R_WH, p_WH, G, b_G, 15);
% stabilityMarginOptimization(kFrictionE, kFrictionH, ...
%         CP_W_e, CN_W_e, CP_H_h, CN_H_h, R_WH, p_WH, e_mode, h_mode);

% parameter sweep
