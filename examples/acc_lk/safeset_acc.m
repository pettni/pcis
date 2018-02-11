mptopt('lpsolver', 'mosek');

% define constants
con = constants;

% Model x^+ = (A + Ap p)x + B u + F + Fp p

% Nominal system
A = [con.f1bar/con.m 0; -1 0];
B = [1/con.m; 0];
F = [-con.f0bar/con.m; con.vl];
XU = Polyhedron('H', [0 0 1 con.Fw_max; 0 0 -1 -con.Fw_min]);

% Parameter dependence
Ap = [0 0; 0 0];
Fp = [-1; 0];

% Parameter bounds
P = Polyhedron('V', [con.nu_max*con.r_max; -con.nu_max*con.r_max]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Compute systems in convex hull %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A_d = eye(2) + con.dt * A;
B_d = con.dt*B;
F_d = con.dt*F;

Ap_d = con.dt * Ap;
Fp_d = con.dt * Fp;

dyn = Dyn(A_d, F_d, B_d, XU, {Ap_d}, {Fp_d}, P);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Compute invariant set %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = Polyhedron('H', [1 0 con.u_max; -1 0 -con.u_min; 0 -1 -con.h_min]);

Cinv = dyn.win_always(X, con.rho_acc, 0, 1);

if ~Cinv.isEmptySet
  poly_A = Cinv.A;
  poly_b = Cinv.b;
  save('results/acc_cinv', 'poly_A', 'poly_b', 'con')
end 
