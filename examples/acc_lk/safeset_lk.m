mptopt('lpsolver', 'gurobi');

% define constants
con = constants;

% Model (z forward speed)
% \dot x = (A + Ap1/z + Ap2*z ) x + B u + E v

A = [0 1 0 0;
     0 0 0 0;
     0 0 0 1;
     0 0 0 0];
B = [0; con.Caf/con.m; 0; con.a*(con.Caf/con.Iz)];
Ev = [0; 0; -1; 0];

XU = Polyhedron('H', [0 0 0 0  1 con.df_max; 
                      0 0 0 0 -1 con.df_max]);

% Parameter dependence
c22 = -((con.Caf+con.Car)/con.m);
c24 = (con.b*con.Car - con.a*con.Caf)/con.m;
c42 = ((con.b*con.Car-con.a*con.Caf)/con.Iz);
c44 = -((con.a^2*con.Caf+con.b^2*con.Car)/con.Iz);

Ap1 = [0 0 0 0;
       0 c22 0 c24;
       0 0 0 0;
       0 c42 0 c44];

Ap2 = [0 0 1  0;
       0 0 0 -1;
       0 0 0  0;
       0 0 0  0];

Fp1 = zeros(4,1);
Fp2 = zeros(4,1);

% Over-approximation of parameters
syms z;
P = compute_hull([1/z, z], [con.u_min con.u_max], 1);

% External disturbance
if con.rd_max > con.rd_min
  V_max = con.rd_min - con.u_min*(con.rd_max-con.rd_min)/(con.u_max-con.u_min);
  V_max_p2 = (con.rd_max-con.rd_min)/(con.u_max-con.u_min);
  XV_V = {[0 0 0 0 0 V_max_p2 V_max], [0 0 0 0 0 -V_max_p2 -V_max]};
else
  XV_V = {[0 0 0 0 0 0 con.rd_max], [0 0 0 0 0 0 -con.rd_max]};
end

% Time discretization
A_d  = eye(4) + con.dt * A;
B_d  = con.dt * B;
Ev_d = con.dt * Ev;

Ap1_d = con.dt*Ap1;
Ap2_d = con.dt*Ap2;

% Dyn object
dyn = Dyn(A_d, [], B_d, XU, ...
          {Ap1_d, Ap2_d}, {Fp1, Fp2}, P, ...
          [], [], [], ...
          Ev_d, XV_V);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Compute invariant set %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

X = Polyhedron('A', [eye(4); -eye(4)], ...
			         'b', [con.y_max; con.nu_max; con.psi_max; con.r_max; ...
                     con.y_max; con.nu_max; con.psi_max; con.r_max]);

Cinv = dyn.win_always(X, con.rho_lk, 0, 1);

if ~Cinv.isEmptySet
  poly_A = Cinv.A;
  poly_b = Cinv.b;
  save('results/lk_cinv', 'poly_A', 'poly_b', 'con')
end 
