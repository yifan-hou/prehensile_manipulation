function [n_av, n_af, R_a, R_a_inv, w_av, Cv, b_C] = hybridServoing(N_all, Nue, G, b_G)

dims.Actualized = 3;
dims.UnActualized = 3;
dims.SlidingFriction = 0;

[n_av, n_af, R_a, w_av] = solvehfvc(dims, N_all, G, ...
    b_G, [], [], [], [], [], 3, false);


if isempty(R_a) || any(isnan(w_av))
    disp('[Hybrid servoing] solvehfvc returns no solution.')
    n_av = []; n_af = [];R_a = [];R_a_inv=[];w_av=[];Cv=[]; b_C=[];
    return;
end

fprintf('[Hybrid Servoing] force dimension: %d\n', n_af);
fprintf('[Hybrid Servoing] velocity dimension: %d\n', n_av);

assert(n_af < 3);
assert(n_af > 0);

% make sure all velocity commands >= 0
R_id_flip = find(w_av < 0);
R_a(n_af + R_id_flip, :) = - R_a(n_af + R_id_flip, :);
w_av(R_id_flip) = - w_av(R_id_flip);

Cv = [zeros(n_av, 3), R_a(n_af + 1: end, :)];
b_C = w_av;
R_a_inv = R_a^-1;
