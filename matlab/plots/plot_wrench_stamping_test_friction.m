clear;clc;

% read friction testing results
fileID = fopen('~/Git/prehensile_manipulation/data/hyper2/test1_friction.txt','r');
A = fscanf(fileID,'%f');

N = A(1);
frictionE_nominal = A(2);
frictionH_nominal = A(3);

frictionE_min = A(4);
frictionE_max = A(5);
frictionH_min = A(6);
frictionH_max = A(7);

X_f = zeros(N*N, 1);
X_i = zeros(N*N, 1);
Y_f = zeros(N*N, 1);
Y_i = zeros(N*N, 1);
f_count = 0;
i_count = 0;

for i = 1:N
    ratio1 = (i-1)/N;
    frictionE = frictionE_min*ratio1 + frictionE_max*(1-ratio1);
    for j = 1:N
        ratio2 = (j-1)/N;
        frictionH = frictionH_min*ratio2 + frictionH_max*(1-ratio2);
        if A(7+(i-1)*N + j) > 0.5
            % feasible
            f_count = f_count+1;
            X_f(f_count) = frictionE;
            Y_f(f_count) = frictionH;
        else
            % infeasible
            i_count = i_count + 1;
            X_i(i_count) = frictionE;
            Y_i(i_count) = frictionH;
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
plot(frictionE_nominal, frictionH_nominal, '.g', 'markersize', 35);
% axis([frictionE_min-0.01, frictionE_max+0.01, frictionH_min-0.01, frictionH_max+0.01]);