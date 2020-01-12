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
kFrictionE = 0.6;

% list of contact points and contact normals
p_e1 = [0.5; 0; 0];
p_e2 = [-0.5; 0; 0];
n_e1 = [0; 1; 0];
n_e2 = [0; 1; 0];

p_O_h1 = [0; 1; 0];
n_O_h1 = [0; -1; 0];

CP_e = [p_e1, p_e2];
CN_e = [n_e1, n_e2];
CP_O_h = [p_O_h1];
CN_O_h = [n_O_h1];

R_WO = eye(3);
p_WO = zeros(3, 1);
%%
%% Goal
%% 0, 1 or 2

% % rrf
% G = [1 0 0 0 0 0];
% b_G = [0.1];
% e_mode = int8([2; 2]); % rr
% h_mode = int8(1); % f

% fsf
G = [1 1 0 0 0 0];
b_G = 0.1;
e_mode = int8([1; 0]); % fs
h_mode = int8(1); % f

% fsf
G = [0 0 1 0 0 0];
b_G = -0.1;
e_mode = int8([1; 0]); % fs
h_mode = int8(1); % f

% wrenchSpaceAnalysis_modeSelection(kFrictionE, kFrictionH, CP_e, CN_e, CP_O_h, CN_O_h, R_WO, p_WO, G, b_G);

wrenchSpaceAnalysis(kFrictionE, kFrictionH, CP_e, CN_e, CP_O_h, CN_O_h, R_WO, p_WO, G, b_G, e_mode, h_mode);
