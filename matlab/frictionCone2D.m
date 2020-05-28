% Compute 2D friction cone edges.
%
% CP: 3xN a set of contact point positions. Last row is all zero.
% CN: 3xN a set of contact normals pointing towards the object. Last row is all zero.
% mu: the coefficient of (static) friction.
% CPF: 3x(2N) a set of contact point positions. Two columns per contact point.
% CNF: 3x(2N) a set of friction cone edges, pointing towards the object. Left edge 1, right edge 1, left edge 2, right edge 2 ...
%
function [CPF, CNF] = frictionCone2D(CP, CN, mu)
N = size(CP,2);
%% CPF
CPF = zeros(3,2*N);
for i=1:N
    CPF(:,(i-1)*2+1:i*2) = CP(:,i)*ones(1,2);
end
%% CNF
CNF = zeros(3,2*N);

% rotate
theta = atan(mu);
c = cos(theta);
s = sin(theta);
RL = [c -s 0; s c 0; 0 0 1]; % positive rotation
RR = [c s 0; -s c 0; 0 0 1];
for i = 1:N
    nn = CN(:, i)/norm(CN(:, i));
    assert(abs(nn(3)) < 1e-7);
    CNF(:,2*i-1) = RL*nn;
    CNF(:,2*i) = RR*nn;
end
assert(norm(CNF(3, :)) < 1e-7);