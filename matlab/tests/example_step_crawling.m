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

kH = 0.1; % object height

xrange = [-0.2 0.2];

% list of contact points and contact normals
n_W_e1 = [0; 1];
n_W_e2 = [0; 1];

p_H_h1 = [0.03; 0];
n_H_h1 = [0; -1];

CN_W_e = [n_W_e1, n_W_e2];
CP_H_h = [p_H_h1];
CN_H_h = [n_H_h1];

R_WH = eye(2);
p_WH = [0; kH];


e_mode = int8([3; 3]); % ll
h_mode = int8(1); % f


% parameter sweep

N = 1000;
px_record = zeros(2, N);
margin_record = zeros(1, N);

for i = 1:N
    k = rand();
    px1 = k*xrange(1) + (1-k)*xrange(2);
    k = rand();
    px2 = k*xrange(1) + (1-k)*px1;

    p_W_e1 = [px1; 0];
    p_W_e2 = [px2; 0];
    CP_W_e = [p_W_e1, p_W_e2];

    margin = computeStabilityMargin(...
            kFrictionE, kFrictionH, CP_W_e, CN_W_e,...
            CP_H_h, CN_H_h, R_WH, p_WH, e_mode, h_mode);

    px_record(:, i) = [px1; px2];
    margin_record(i) = margin;
end

% plot
id_stable = margin_record > 1e-7;


cherry=[153/255,15/255,2/255];
fern=[92/255,188/255,99/255];
onyx=[3/255,1/255,6/255];


figure(1);clf(1); hold on;
plot(px_record(1,id_stable), px_record(2,id_stable), '.', 'color', fern);
plot(px_record(1,~id_stable), px_record(2,~id_stable), '.', 'color', cherry);
plot(0.05, -0.05, '.', 'markersize', 15, 'color', onyx);

[hLegend, icons] = legend('Stable', 'Unstable', 'box','off');
hxl = xlabel('Right Corner Location (m)');
hyl = ylabel('Left Corner Location (m)');

icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
set(icons,'MarkerSize',20);

set(gca, 'FontName', 'Times New Roman')
set([hLegend, hxl, hyl], 'FontName', 'Times New Roman');
set(hLegend,'FontSize',10);
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1);

set(gcf, 'Position',  [200, 300, 400, 250])
