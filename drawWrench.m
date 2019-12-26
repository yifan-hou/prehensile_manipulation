% Draw 3D wrenches

function drawWrench(W, color)

origin = [0 0 0];
data = [origin; W'];
if size(W, 2) >= 3
    K = convhull(data(:,1), data(:,2), data(:,3));
    h = trisurf(K, data(:,1), data(:,2), data(:,3), 'Facecolor',color);
else
    h = patch(data(:,1), data(:,2), data(:,3), color);
end

h.FaceColor = color;
h.FaceAlpha = 0.5;
% axis equal;


xlabel('X');
ylabel('Y');
zlabel('Z');

end