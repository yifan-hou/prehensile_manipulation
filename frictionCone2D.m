% 16-741 Mechanics of Manipulation, Fall 2015
% Author: Sung Kyun Kim (kimsk@cs.cmu.edu)
%
% N: the number of contact points; scalar
% CP: a set of contact point positions [[pix; piy; piz] ...]; 3xN matrix
% CN: a set of inward-pointing directions of contact normals [[nix; niy; niz] ...]; 3xN matrix
% mu: the coefficient of (static) friction; scalar
% M: the number of side facets of a linearized polyhedral friction cone; scalar
% CPF: a set of contact point positions of edges of polyhedral friction cones [[pijx; pijy; pijz] ...]; 3x(NM) matrix
% CNF: a set of inward-pointing directions of edges of polyhedral friction cones [[sijx; sijy; sijz] ...]; 3x(NM) matrix
%
% Examples:
% i;                            % index for i-th contact normal
% j;                            % index for j-th edge of polyhedral friction cone of i-th contact normal
% CPi = CP(1:3, i);             % i-th contact point position
% CNi = CN(1:3, i);             % i-th contact normal direction
% CPFij = CPF(1:3, M*(i-1)+j);  % contact point position of j-th edge of polyhedral friction cone of i-th contact normal
% CNFij = CNF(1:3, M*(i-1)+j);  % contact normal direction of j-th edge of polyhedral friction cone of i-th contact normal

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
RL = [c -s 0; s c 0; 0 0 1];
RR = [c s 0; -s c 0; 0 0 1];
for i = 1:N
    nn = CN(:, i)/norm(CN(:, i));
    assert(abs(nn(3)) < 1e-7);
    CNF(:,2*i-1) = RL*nn;
    CNF(:,2*i) = RL*nn;
end
assert(norm(CNF(3, :)) < 1e-7);