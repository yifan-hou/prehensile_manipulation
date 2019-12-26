% Double Description Method
% Given a homogeneous polyhedral convex cone described by {x: Ax <= 0},
% Find its ray description {y: y=Rx, x >= 0}
% Note: this representation may not be the minimal representation. Some columns of R might be redundant.
%
%   Examples:
%       A = [0 1 0; 1 0 0];
%       R = DoubleDescription(A);
%
% Implemented based on
%   https://inf.ethz.ch/personal/fukudak/lect/pclect/notes2014/PolyComp2014.pdf
%   Lemma 9.1
function [R] = DoubleDescription(A)

[kRowsA, kDim] = size(A);
% Trivial initial DD pair
R = [eye(kDim) -ones(kDim, 1)]; % initially, R expands the whole space

% ID_ALL = 1:kRowsA;
% while length(K) < kRowsA
for i = 1:kRowsA
    % select any index i from {1, 2, ..., kRowsA}\K
    Ai = A(i, :);
    % this new hyperplane, Ai, divides rays in R into three groups
    projection = Ai*R;
    J_minus = projection < -1e-10;
    J_plus = projection > 1e-10;
    J_zero = ~(J_minus | J_plus);
    % new R is composed of rays of J_minus, rays of J_zero, and some new vectors
    R1 = R(:, J_minus | J_zero);
    R2 = [];
    id_plus = find(J_plus);
    id_minus = find(J_minus);
    for j = 1:sum(J_plus)
        for j_ = 1:sum(J_minus)
            r_jj_ = (Ai*R(:, id_plus(j)))*R(:, id_minus(j_)) - (Ai*R(:, id_minus(j_)))*R(:, id_plus(j));
            R2 = [R2 r_jj_];
        end
    end
    % update
    R = [R1 R2];
end
