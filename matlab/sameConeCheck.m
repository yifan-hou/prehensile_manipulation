% cone1: 3xn1
% cone2: 3xn2
function is_same = sameConeCheck(cone1, cone2)
    is_same = false;
    n1 = size(cone1, 2);
    n2 = size(cone2, 2);
    cone1 = normc(cone1);
    cone2 = normc(cone2);
    if n1 ~= n2
        return;
    end

    TOL = 1e-7;
    orders = perms(n1:-1:1);
    for i = 1:size(orders, 1)
        if norm(cone1 - cone2(:, orders(i, :))) < TOL
            is_same = true;
            return;
        end
    end
end