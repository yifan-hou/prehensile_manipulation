function N = getJacobianFromContacts(mode_e, mode_h, Jac_e, Jac_h)

Ne = length(mode_e);
Nh = length(mode_h);

% 0:separation 1:fixed 2/3: sliding
nRows = 2*sum(mode_e == 1) + sum(mode_e >= 2) + ...
        2*sum(mode_h == 1) + sum(mode_h >= 2);

N = zeros(nRows, 6);

Ncount = 0;
for i = 1:Ne
    if mode_e(i) == 0
        continue;
    end
    Ncount = Ncount + 1;
    N(Ncount, :) = Jac_e(2*i-1, :);
    if mode_e(i) == 1
        Ncount = Ncount + 1;
        N(Ncount, :) = Jac_e(2*i, :);
    end
end

for i = 1:Nh
    if mode_h(i) == 0
        continue;
    end
    Ncount = Ncount + 1;
    N(Ncount, :) = Jac_h(2*i-1, :);
    if mode_h(i) == 1
        Ncount = Ncount + 1;
        N(Ncount, :) = Jac_h(2*i, :);
    end
end