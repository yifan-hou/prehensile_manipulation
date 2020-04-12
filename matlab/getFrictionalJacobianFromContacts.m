% N: bilateral constraints. N*v = 0. Includes:
%   1. vn = 0 for sliding
%   2. [vn, vt] = 0 for sticking
% Nu: unilateral constraints. Nu*v > 0. Includes:
%   1. vn >=0 for separation.
%   2. vt >= 0 for left sliding, vt <= 0 for right sliding.
%
% id_e: id of rows of Jac_e that goes to Je_
% id_h: id of rows of Jac_h that goes to Jh_
function [Je_, Jh_, id_e, id_h] = getFrictionalJacobianFromContacts(mode_e, mode_h, Jac_e, Jac_h)

Ne = length(mode_e);
Nh = length(mode_h);
% 0:separation 1:fixed 2/3: sliding
nRows_e = 2*sum(mode_e == 1) + sum(mode_e >= 2);
nRows_h = 2*sum(mode_h == 1) + sum(mode_h >= 2);

Je_ = zeros(nRows_e, 3);
Jh_ = zeros(nRows_h, 3);
id_e = zeros(nRows_e, 1);
id_h = zeros(nRows_h, 1);

e_count = 0;
for i = 1:Ne
    if mode_e(i) == 1
        e_count = e_count + 2;
        Je_(e_count-1:e_count, :) = Jac_e(2*i - 1:2*i, :);
        id_e(e_count-1:e_count) = 2*i - 1:2*i;
    elseif mode_e(i) == 2
        % right sliding, left edge of friction cone
        e_count = e_count + 1;
        Je_(e_count, :) = Jac_e(2*i-1, :);
        id_e(e_count) = 2*i-1;
    elseif mode_e(i) == 3
        % left sliding, right edge of friction cone
        e_count = e_count + 1;
        Je_(e_count, :) = Jac_e(2*i, :);
        id_e(e_count) = 2*i;
    end
end

h_count = 0;
for i = 1:Nh
    if mode_h(i) == 1
        h_count = h_count + 2;
        Jh_(h_count-1:h_count, :) = Jac_h(2*i - 1:2*i, :);
        id_h(h_count-1:h_count) = 2*i - 1:2*i;
    elseif mode_h(i) == 2
        % right sliding, left edge of friction cone
        h_count = h_count + 1;
        Jh_(h_count, :) = Jac_h(2*i-1, :);
        id_h(h_count) = 2*i-1;
    elseif mode_h(i) == 3
        % left sliding, right edge of friction cone
        h_count = h_count + 1;
        Jh_(h_count, :) = Jac_h(2*i, :);
        id_h(h_count) = 2*i;
    end
end