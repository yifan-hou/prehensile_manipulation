% 16-741 Mechanics of Manipulation, Fall 2015
% Author: Sung Kyun Kim (kimsk@cs.cmu.edu)
%
% This function draws given contact normals in a 3D plot
% Forces are drawn at contact points, moments are drawn at the origin
%
% N: the number of contact points; scalar
% CP: a set of contact point positions [[pix; piy; piz] ...]; 3x(NM) matrix
% W: a set of normalized contact screws [[cix; ciy; ciz; c0ix; c0iy; c0iz] ...]; 6x(NM) matrix
% M: the number of side facets of a linearized polyhedral friction cone; scalar
%
% Examples:
% drawContactScrew(CP, W);      % for frictionless point contact normals
% drawContactScrew(CP, W, M);   % for frictional point contact normals

function drawContactScrew(CP, W, M)

% check input arguments
if nargin == 2
    M = 1;
end
if size(CP,2) ~= size(W,2)
    error('The number of columns of CP and W should be same!');
end
if mod(size(W,2), M) ~= 0
    error('The number of column vectors of W should be a multiple of M!');
end
N = size(W,2)/M;    % the number of contact points

% plot
hold on; grid on;
O = zeros(1,M);
for i = 1:N
    CPi = CP(:, M*(i-1)+1:M*i);
    Wi = W(:, M*(i-1)+1:M*i);
    coli = rand(1,3);   % color for each friction cone

    % draw contact points
    plot3(CPi(1),CPi(2),CPi(3), 'kp','markersize',12);
    plot3(O,O,O, 'ko','markersize',10);

    % draw contact normals (the point of action at the end of arrow)
    quiver3(CPi(1,:),CPi(2,:),CPi(3,:), Wi(1,:),Wi(2,:),Wi(3,:), 'autoscale','off','linewidth',2,'color',coli);     % force applied at contact point
    quiver3(O,O,O, Wi(4,:),Wi(5,:),Wi(6,:), 'autoscale','off','linewidth',2,'color',coli);                          % moments about the origin
end
title('Contact normals in 3D space');
xlabel('x');
ylabel('y');
zlabel('z');
axis equal;
view([10 20]);
