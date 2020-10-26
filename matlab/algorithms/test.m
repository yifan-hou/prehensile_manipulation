J = [  0.106067,   0.106067,    0.02175,          0,          0,          0;
   0.106067,   0.106067, -0.0217499,          0,          0,          0;
  -0.106066,  -0.106066,         -0,   0.106067,   0.106067,          0;
   0.106067,  -0.106066,         -0,  -0.106066,   0.106067,          0];


n_a = 3;
n_u = 3;
n   = 6;

M = [zeros(n_u, n_a); eye(n_a)];
U = null(J)';

if isempty(U)
    % no null space. The environment already fully constrained the problem
    return;
end

U_bar = U*M;
U_hat = chol(U_bar'*U_bar + 1e-15*eye(n_a));

TOL = 0.0001; % this tolerance has to be big.
% for i = 1:n_a
%     if (norm(U_hat(n_a-i + 1,:)) > TOL)
%         U_hat = U_hat(1:n_a-i+1, :);
%         break;
%     end
% end
id = [];
for i = 1:n_a
    if (norm(U_hat(n_a-i + 1,:)) > TOL)
        id = [id n_a-i+1];
    end
end
U_hat = U_hat(id, :);


n_av = size(U_hat,1);
n_af = n_a - n_av;

% C_bar = (U_hat\eye(n_av))';
C_bar = lsqminnorm(U_hat, eye(n_av))';
R_a = [null(C_bar)';
        C_bar];
T = blkdiag(eye(n_u), R_a);
C = [zeros(n_av, n_u) C_bar];
% solve for b_C

JC = [J; C];