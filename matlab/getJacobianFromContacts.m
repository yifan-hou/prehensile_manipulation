% N: bilateral constraints. N*v = 0. Includes:
%   1. vn = 0 for sliding
%   2. [vn, vt] = 0 for sticking
% Nu: unilateral constraints. Nu*v > 0. Includes:
%   1. vn >=0 for separation.
%   2. vt >= 0 for left sliding, vt <= 0 for right sliding.
% Nue: the equality constraints that are actually active unilateral constraints
%   1. vn >= 0 for sticking and sliding
function [N, Nu, Nue] = getJacobianFromContacts(mode_e, mode_h, Jac_e, Jac_h)

Ne = length(mode_e);
Nh = length(mode_h);

% 0:separation 1:fixed 2/3: sliding
nRows = 2*sum(mode_e == 1) + sum(mode_e >= 2) + ...
        2*sum(mode_h == 1) + sum(mode_h >= 2);
nRows_u = sum(mode_e == 0) + sum(mode_e >= 2) + ...
          sum(mode_h == 0) + sum(mode_h >= 2);
nRows_ue = sum(mode_e >= 1) + sum(mode_h >= 1);

N = zeros(nRows, 6);
Nu = zeros(nRows_u, 6);
Nue = zeros(nRows_ue, 6);

Ncount = 0;
Ncount_u = 0;
Ncount_ue = 0;
for i = 1:Ne
    if mode_e(i) == 0
        % separating
        Ncount_u = Ncount_u + 1;
        Nu(Ncount_u, :) = Jac_e(2*i-1, :);
        continue;
    end
    Ncount = Ncount + 1;
    N(Ncount, :) = Jac_e(2*i-1, :);
    Ncount_ue = Ncount_ue + 1;
    Nue(Ncount_ue, :) = Jac_e(2*i-1, :);
    if mode_e(i) == 1
        % fixed
        Ncount = Ncount + 1;
        N(Ncount, :) = Jac_e(2*i, :);
    elseif mode_e(i) == 2
        % right sliding
        Ncount_u = Ncount_u + 1;
        Nu(Ncount_u, :) = -Jac_e(2*i, :);
    else
        % left sliding
        Ncount_u = Ncount_u + 1;
        Nu(Ncount_u, :) = Jac_e(2*i, :);
    end
end

for i = 1:Nh
    if mode_h(i) == 0
        % separating
        Ncount_u = Ncount_u + 1;
        Nu(Ncount_u, :) = Jac_h(2*i-1, :);
        continue;
    end
    Ncount = Ncount + 1;
    N(Ncount, :) = Jac_h(2*i-1, :);
    Ncount_ue = Ncount_ue + 1;
    Nue(Ncount_ue, :) = Jac_h(2*i-1, :);
    if mode_h(i) == 1
        % fixed
        Ncount = Ncount + 1;
        N(Ncount, :) = Jac_h(2*i, :);
    elseif mode_h(i) == 2
        % right sliding
        Ncount_u = Ncount_u + 1;
        Nu(Ncount_u, :) = -Jac_h(2*i, :);
    else
        % left sliding
        Ncount_u = Ncount_u + 1;
        Nu(Ncount_u, :) = Jac_h(2*i, :);
    end
end
