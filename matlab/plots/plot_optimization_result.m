clear;clc;
load generated/llf.mat

m_llf = margin_record;
p_llf = p_record;

load generated/sff.mat

m_sff = margin_record;
p_sff = p_record;

% figure(1);clf(1);hold on;
% subplot(2,1,1); hold on;
% hline1 = plot(m_sff,'-','Color', [0 0 .5]);
% hline2 = plot(m_llf,'-','Color', [0 .5 .5]);
% hLegend = legend('Pivoting','Sliding', 'box','off');
% hTitle = title('Stability Margin (Unitless)');
%
% set(gca, 'FontName', 'Times New Roman')
% set([hTitle], 'FontName', 'Times New Roman')
% set([hline1 hline2], 'LineWidth', 1.5)
% set([hLegend, gca], 'FontSize', 8)
% set(hTitle, 'FontSize', 10, 'FontWeight' , 'bold')
% % Adjust axes properties
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
%     'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:0.2:1, ...
%     'LineWidth', 1)
%
% subplot(2,1,2); hold on;
% hline1 = plot(p_sff,'-','Color', [0 0 .5]);
% hline2 = plot(p_llf,'-','Color', [0 .5 .5]);
% hLegend = legend('Pivoting','Sliding', 'box','off');
% hTitle = title('Finger Position X (m)');
% hxl = xlabel('Time steps');
%
% axis([0 200 -0.05 0.03]);
% set(gca, 'FontName', 'Times New Roman')
% set([hTitle], 'FontName', 'Times New Roman')
% set([hline1 hline2], 'LineWidth', 1.5)
% set([hLegend, gca], 'FontSize', 8)
% set(hTitle, 'FontSize', 10, 'FontWeight' , 'bold')
% % Adjust axes properties
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
%     'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', -0.05:0.02:0.03, ...
%     'LineWidth', 1)
%
% set(gcf, 'Position',  [200, 200, 250, 350])

figure(1);clf(1);hold on;
subplot(1,2,1); hold on;
hline1 = plot(m_sff,'-','Color', [0 0 .5]);
hline2 = plot(m_llf,'-','Color', [0 .5 .5]);
hLegend = legend('Pivoting','Sliding', 'box','off');
hTitle = title('Stability Margin (Unitless)');
hx1 = xlabel('Time steps');

set(gca, 'FontName', 'Times New Roman')
set([hTitle], 'FontName', 'Times New Roman')
set([hline1 hline2], 'LineWidth', 1.5)
set([hLegend, gca], 'FontSize', 8)
set(hTitle, 'FontSize', 10, 'FontWeight' , 'bold')
% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', 0:0.2:1, ...
    'LineWidth', 1)

subplot(1,2,2); hold on;
hline1 = plot(p_sff,'-','Color', [0 0 .5]);
hline2 = plot(p_llf,'-','Color', [0 .5 .5]);
hLegend = legend('Pivoting','Sliding', 'box','off');
hTitle = title('Finger Position X (m)');
hx2 = xlabel('Time steps');

axis([0 200 -0.05 0.03]);
set(gca, 'FontName', 'Times New Roman')
set([hTitle], 'FontName', 'Times New Roman')
set([hline1 hline2], 'LineWidth', 1.5)
set([hLegend, gca], 'FontSize', 8)
set(hTitle, 'FontSize', 10, 'FontWeight' , 'bold')
% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'YTick', -0.05:0.02:0.03, ...
    'LineWidth', 1)

% set(gcf, 'Position',  [200, 200, 250, 350])
