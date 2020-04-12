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

function [CPF, CNF] = frictionCone(CP, CN, mu, M)
N = size(CP,2);
%% CPF
CPF = zeros(3,N*M);
for i=1:N
    CPF(:,(i-1)*M+1:i*M) = CP(:,i)*ones(1,M);
end
%% CNF
CNF = zeros(3,N*M);
% module, a cone with axis [0 0 1]'
module = ones(3,M);
for i = 1:M
    theta = (i-1)/M * 2*pi;
    module(1,i) = mu*cos(theta);
    module(2,i) = mu*sin(theta);
end
% rotate
for i = 1:N
    nn = CN(:,i)/norm(CN(:,i));
    axis = cross([0 0 1]',nn);
    if norm(axis) < 1e-7
        axis = [1 0 0]';
    end
    axis = axis/norm(axis);
    Ncross = [0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];
    R = eye(3) + norm(nn(1:2))*Ncross + (1 - nn(3))*Ncross*Ncross;
    CNF(:,(i-1)*M+1:i*M) = R*module;
end