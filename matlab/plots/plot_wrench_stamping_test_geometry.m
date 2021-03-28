clear;clc;

% read friction testing results
fileID = fopen('~/Git/prehensile_manipulation/data/test2_geometrical.txt','r');
A = fscanf(fileID,'%f');

N = A(1);
L1_nominal = A(2);
L2_nominal = A(3);

L1_min = A(4);
L1_max = A(5);
L2_min = A(6);
L2_max = A(7);

X_f = zeros(N*N, 1);
X_i = zeros(N*N, 1);
Y_f = zeros(N*N, 1);
Y_i = zeros(N*N, 1);
f_count = 0;
i_count = 0;

for i = 1:N
    ratio1 = (i-1)/N;
    L1 = L1_min*ratio1 + L1_max*(1-ratio1);
    for j = 1:N
        ratio2 = (j-1)/N;
        L2 = L2_min*ratio2 + L2_max*(1-ratio2);
        if A(7+(i-1)*N + j) > 0.5
            % feasible
            f_count = f_count+1;
            X_f(f_count) = L1;
            Y_f(f_count) = L2;
        else
            % infeasible
            i_count = i_count + 1;
            X_i(i_count) = L1;
            Y_i(i_count) = L2;
        end
    end
end

X_i = X_i(1:i_count);
Y_i = Y_i(1:i_count);
X_f = X_f(1:f_count);
Y_f = Y_f(1:f_count);

figure(1);clf(1);hold on;
plot(X_f, Y_f, '.b', 'markersize', 10);
plot(X_i, Y_i, '.r', 'markersize', 10);
plot(L1_nominal, L2_nominal, '.g', 'markersize', 35);
% axis([L1_min-0.01, L1_max+0.01, L2_min-0.01, L2_max+0.01]);