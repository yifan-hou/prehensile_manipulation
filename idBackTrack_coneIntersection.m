% Ca: cone a, 3xNa
% Cb: cone b, 3xNb
% Ci: cone intersection, 3xNi
% idm: logical matrix
%       a1 a2 a3 ... b1 b2 ..
%   i1: 1  0  0      0  0    (vertex of Ca)
%   i2: 1  1  0      1  1    (Intersection of edges)
function idm = idBackTrack_coneIntersection(Ca, Cb, Ci)

Na = size(Ca, 2);
Nb = size(Cb, 2);
Ni = size(Ci, 2);

Ca = normc(Ca);
Cb = normc(Cb);
Ci = normc(Ci);

TOL = 1e-8;

idm = false(Ni, Na+Nb);
for idi = 1:Ni
    % vertex check
    ray_i = Ci(:, idi);
    found = false;
    % cone a
    for ida = 1:Na
        if 1 - (ray_i'*Ca(:, ida)) < TOL
            found = true;
            break;
        end
    end
    if found
        idm(idi, ida) = true;
        continue;
    end
    % cone b
    for idb = 1:Nb
        if 1 - (ray_i'*Cb(:, idb)) < TOL
            found = true;
            break;
        end
    end
    if found
        idm(idi, Na+idb) = true;
        continue;
    end

    % edge intersection check
    for ida1=1:Na
        height = cross(Ca(:, ida1), ray_i);
        for ida2 = ida1+1:Na
            if abs(height'*Ca(:, ida2)) < TOL
                found = true;
                break;
            end
        end
        if found
            break;
        end
    end
    if ~found
        error('[idBackTrack_coneIntersection] error: intersection not found');
    end
    found = false;
    for idb1=1:Nb
        height = cross(Cb(:, idb1), ray_i);
        for idb2 = idb1+1:Nb
            if abs(height'*Cb(:, idb2)) < TOL
                found = true;
                break;
            end
        end
        if found
            break;
        end
    end
    if ~found
        error('[idBackTrack_coneIntersection] error: intersection not found');
    end
    idm(idi, [ida1 ida2 Na+idb1 Na+idb2]) = true;
end
