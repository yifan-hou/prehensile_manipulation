%!
%! Calculates the stability margin for cone intersection. cone_inside is inside
%! of cone_outside. The function computes how deep cone_inside is inside of
%! cone_outside. Formula:
%!
%! margin_ = min_j {max_i Distance(Ej, ai)} where Ej is the jth edge of
%! cone_outside, ai is the ith ray of cone_inside.
%!
%! cone_outside:  3 x n1, the n1 rays of the outside cone.
%! cone_inside:   3 x n2, the n2 rays of the inside cone.
%!
%! margin_:     The stability margin. The more stable, the larger margin_ should
%!             be. When margin_ is zero, the two cones only intersect on the
%!             boundary. When margin_ is negative, the two cones only intersect
%!             at the origin.
%!
%! i:   index of the ray ai in cone_outside.
%!
%! j, k:    index of rays in cone_outside that forms edge Ej.
function [margin_, i, j, k] = computeStabilityMargin(cone_outside, cone_inside)

[d, n1] = size(cone_outside);
n2 = size(cone_inside, 2);
assert(d == 3);

cone_outside = normc(cone_outside);
cone_inside = normc(cone_inside);

origin = [0 0 0]';
if n1 >= 3
    K = convhull([origin cone_outside]');
    tri = K(any(K==1, 2), :);
    assert(size(tri, 1) == n1);
    center = mean(cone_outside, 2);

    i = 0;
    j = 0;
    k = 0;
    margin_ = inf;
    for ej = 1:n1
        ids = tri(ej, :) - 1;
        ids(ids == 0) = [];
        assert(length(ids) == 2);
        facet = cone_outside(:, ids);
        n_in = cross(facet(:, 1), facet(:, 2));
        n_in = n_in/norm(n_in);
        if n_in'*center < 0
            n_in = -n_in;
        end
        dists = n_in'*cone_inside;
        [new_dist, new_i] = max(dists);
        if new_dist < margin_
            margin_ = new_dist;
            i = new_i;
            j = ids(1);
            k = ids(2);
        end
    end
elseif n1 == 2
    j = 1;
    k = 2;
    z = cross(cone_outside(:, 1), cone_outside(:, 2));
    z = z/norm(z);

    angs1 = zeros(1, n2);
    angs2 = zeros(1, n2);
    for ei = 1:n2
        angs1(ei) = abs(angBTVec(cone_outside(:, 1), cone_inside(:, ei), z));
        angs2(ei) = abs(angBTVec(cone_outside(:, 2), cone_inside(:, ei), -z));
    end
    [max1, i1] = max(angs1);
    [max2, i2] = max(angs2);
    if max1 < max2
        margin_ = max1;
        i = i1;
    else
        margin_ = max2;
        i = i2;
    end
else
    % n1 == 1
    % note since cone_inside has to be inside of cone_outside, then cone_inside
    % must be the same ray as cone_outside
    margin_ = 0;
    i = 1;
    j = 1;
    k = [];
end
