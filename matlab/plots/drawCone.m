function drawCone(W, color, addOrigin, texts)

W = normalizeByCol(W);
if addOrigin
    origin = [0 0 0];
    data = [origin; W'];
else
    data = W';
end


% for i = 1:size(W,2)
%     quiver3(0,0,0, W(1, i), W(2, i), W(3, i), 'linewidth',3, 'markersize', 15, 'markerEdgeColor', color);
% end
plot3(W(1, :), W(2, :), W(3, :), '.', 'markersize', 20, 'markerEdgeColor', color);

rw = rank(bsxfun(@minus, data', data(1,:)'));

if rw >= 3
    % 3D polytope
    K = convhull(data(:,1), data(:,2), data(:,3));
    h = trisurf(K, data(:,1), data(:,2), data(:,3));
    h.FaceAlpha = 0.5;
elseif rw == 2
    % 2D plane
    if size(data, 1) >= 3
        [U,S]=svd(   bsxfun(@minus,data,mean(data)),   0);
        K = convhull(U*S(:,1:2));
        h = patch(data(K, 1), data(K, 2), data(K, 3), color);
    else
        h = patch(data(:, 1), data(:, 2), data(:, 3), color);
    end
    h.FaceAlpha = 0.95;
else
    % line or ray
    h = patch(data(:, 1), data(:, 2), data(:, 3), color);
    h.LineWidth = 3;
end

h.FaceColor = color;
h.EdgeColor = color;

% axis equal;

if nargin >= 4
    text(mean(W(1,:)),mean(W(2,:))+0.1,mean(W(3,:)),['   ' ...
        texts], 'HorizontalAlignment','center','FontSize',14,'FontName','Times New Roman');
end

hxl = xlabel('F_X(N)');
hyl = ylabel('F_Y(N)');
hzl = zlabel('T_Z(N)');
grid on;
set([hxl hyl hzl], 'FontName', 'Times New Roman', 'FontWeight','bold')
axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);

set(gca, 'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], 'ZColor', [.3 .3 .3],...
    'LineWidth', 1)
% axis equal

end