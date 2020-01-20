% Draw 3D wrenches

function drawWrench(W, color, addOrigin, texts)

if addOrigin
    origin = [0 0 0];
    data = [origin; W'];
else
    data = W';
end

if rank(W) >= 3
    % 3D polytope
    K = convhull(data(:,1), data(:,2), data(:,3));
    h = trisurf(K, data(:,1), data(:,2), data(:,3));
    h.FaceAlpha = 0.5;
elseif rank(W) == 2
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

plot3(W(1, :), W(2, :), W(3, :), '.', 'markersize', 20, 'markerEdgeColor', color);

if nargin >= 4
    text(mean(W(1,:)),mean(W(2,:)),mean(W(3,:)),['   ' ...
        texts], 'HorizontalAlignment','left','FontSize',18);
end

xlabel('X');
ylabel('Y');
zlabel('Z');

axis([-1.1 1.1 -1.1 1.1 -1.1 1.1]);

axis equal

end