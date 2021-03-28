% N: bilateral constraints. N*v = 0. Includes:
%   1. vn = 0 for sliding
%   2. [vn, vt] = 0 for sticking
% Nu: unilateral constraints. Nu*v > 0. Includes:
%   1. vn >=0 for separation.
%   2. vt >= 0 for left sliding, vt <= 0 for right sliding.
% Nue: the equality constraints that are actually active unilateral constraints
%   1. vn >= 0 for sticking and sliding
function [J, normal_ids] = getJacobianFromContacts3D(mode_e, mode_h, Normal_e, Normal_h, Tangent_e, Tangent_h)

Ne = length(mode_e);
Nh = length(mode_h);
kDim = size(Normal_e, 2);

% 0:separation 1:fixed 2/3: sliding
% 1: sticking   0: sliding
nRowsE = 3*sum(mode_e == 1) + sum(mode_e == 0);
nRowsH = 3*sum(mode_h == 1) + sum(mode_h == 0);

Je = zeros(nRowsE, 6);
Jh = zeros(nRowsH, 6);
NEcount = 1;
NHcount = 1;

normal_ids = [];

for i = 1:Ne
    Je(NEcount, :) = Normal_e(i, :);
    normal_ids = [normal_ids NEcount];
    NEcount = NEcount + 1;
    if mode_e(i) == 1
        Je(NEcount:NEcount + 1, :) = Tangent_e(2*i-1:2*i, :);
        NEcount = NEcount + 2;
    end
end

for i = 1:Nh
    Jh(NHcount, :) = Normal_h(i, :);
    normal_ids = [normal_ids NHcount + nRowsE];
    NHcount = NHcount + 1;
    if mode_h(i) == 1
        Jh(NHcount:NHcount + 1, :) = Tangent_h(2*i-1:2*i, :);
        NHcount = NHcount + 2;
    end
end

J = [Je, zeros(nRowsE, kDim); -Jh, Jh];
