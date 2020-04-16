clc;clear;
% Parameters
kFrictionE = 0.5;

% list of contact points and contact normals
p1 = [0; 1; 0];
p2 = [0; 0; 0];
p3 = [0; 0; 0];
p4 = [1; 0; 0];

n1 = [1; 0; 0];
n2 = [1; 0; 0];
n3 = [0; 1; 0];
n4 = [0; 1; 0];

CP = [p1, p2, p3, p4];
CN = [n1, n2, n3, n4];

% friction cone approximation
[CPF, CNF] = frictionCone2D(CP, CN, kFrictionE);

% contact screws
W = contactScrew(CPF, CNF);

% draw contact screws
figure(1); clf(1); hold on;
drawCone(W([1 2 6], :), 'b');
plot3(W(1, 1), W(2, 1), W(6, 1), 'r*', 'markersize', 15);
plot3(W(1, 2), W(2, 2), W(6, 2), 'g*', 'markersize', 15);
plot3(W(1, 3), W(2, 3), W(6, 3), 'ro', 'markersize', 15);
plot3(W(1, 4), W(2, 4), W(6, 4), 'go', 'markersize', 15);
plot3(W(1, 5), W(2, 5), W(6, 5), 'r.', 'markersize', 15);
plot3(W(1, 6), W(2, 6), W(6, 6), 'g.', 'markersize', 15);
plot3(W(1, 7), W(2, 7), W(6, 7), 'rx', 'markersize', 15);
plot3(W(1, 8), W(2, 8), W(6, 8), 'gx', 'markersize', 15);

view(0, 90);
