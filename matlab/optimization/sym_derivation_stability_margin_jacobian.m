clear;clc;
Ci = sym('Ci', [3,1], 'real');
Cj = sym('Cj', [3,1], 'real');
Ck = sym('Ck', [3,1], 'real');
C1 = sym('C1', [3,1], 'real');
C2 = sym('C2', [3,1], 'real');
C3 = sym('C3', [3,1], 'real');
C4 = sym('C4', [3,1], 'real');

% norm_vec = @(v) sqrt(v'*v);

Cij = cross_vec(Ci, Cj);
C12 = cross_vec(C1, C2);
C34 = cross_vec(C3, C4);
C1234 = cross_vec(C12, C34);

Phi_facet_vertex = asin(Cij'*Ck);
Phi_facet_intersection = asin(Cij'*C1234);

Phi_edge_vertex = acos(Ci'*Ck);
Phi_edge_intersection = acos(Ci'*C1234);


%% ---------------------------------------------------------------
%           calculate derivatives
% ---------------------------------------------------------------
deriv_ijk  = @(f) [
        diff(f,'Ci1'), diff(f,'Ci2'), diff(f,'Ci3'); ...
        diff(f,'Cj1'), diff(f,'Cj2'), diff(f,'Cj3'); ...
        diff(f,'Ck1'), diff(f,'Ck2'), diff(f,'Ck3')];
deriv_ij1234  = @(f) [
        diff(f,'Ci1'), diff(f,'Ci2'), diff(f,'Ci3'); ...
        diff(f,'Cj1'), diff(f,'Cj2'), diff(f,'Cj3'); ...
        diff(f,'C11'), diff(f,'C12'), diff(f,'C13'); ...
        diff(f,'C21'), diff(f,'C22'), diff(f,'C23'); ...
        diff(f,'C31'), diff(f,'C32'), diff(f,'C33'); ...
        diff(f,'C41'), diff(f,'C42'), diff(f,'C43')];
deriv_ik  = @(f) [
        diff(f,'Ci1'), diff(f,'Ci2'), diff(f,'Ci3'); ...
        diff(f,'Ck1'), diff(f,'Ck2'), diff(f,'Ck3')];
deriv_i1234  = @(f) [
        diff(f,'Ci1'), diff(f,'Ci2'), diff(f,'Ci3'); ...
        diff(f,'C11'), diff(f,'C12'), diff(f,'C13'); ...
        diff(f,'C21'), diff(f,'C22'), diff(f,'C23'); ...
        diff(f,'C31'), diff(f,'C32'), diff(f,'C33'); ...
        diff(f,'C41'), diff(f,'C42'), diff(f,'C43')];


jac_phi_facet_vertex = deriv_ijk(Phi_facet_vertex);
jac_phi_facet_intersection = deriv_ij1234(Phi_facet_intersection);
jac_phi_edge_vertex = deriv_ik(Phi_edge_vertex);
jac_phi_edge_intersection = deriv_i1234(Phi_edge_intersection);


jac_phi_ijk = (jac_phi_facet_vertex');
jac_phi_ij1234 = (jac_phi_facet_intersection');
jac_phi_ik = (jac_phi_edge_vertex');
jac_phi_i1234 = (jac_phi_edge_intersection');

disp('Done. Generating file:');

%% ---------------------------------------------------------------
%           write to file
% ---------------------------------------------------------------
matlabFunction(jac_phi_ijk, 'File', 'generated/jac_phi_ijk', 'vars', {Ci, Cj, Ck});
matlabFunction(jac_phi_ij1234, 'File', 'generated/jac_phi_ij1234', 'vars', {Ci, Cj, C1, C2, C3, C4});

matlabFunction(jac_phi_ik, 'File', 'generated/jac_phi_ik', 'vars', {Ci, Ck});
matlabFunction(jac_phi_i1234, 'File', 'generated/jac_phi_i1234', 'vars', {Ci, C1, C2, C3, C4});

disp('All done');
