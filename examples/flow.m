% Finding a controlled invariant set via a one-shot method

% System matrices, state = [x; dx; ddx]
A = eye(5);

F = [0; 0; 0; 0; 4];

B = [-1  0  0  0  1  0;
      1 -1  0  0  0  0;
      0  1 -1  0  0 -1;
      0  0  1 -1  0  0;
      0  0  0  1 -1  0];

nx = 5;
nu = 6;

H_xu = [zeros(nu, nx) -eye(nu) zeros(nu, 1); 
        -1  0  0  0  0    1  0  0  0  0  0      0;
         0 -1  0  0  0    0  1  0  0  0  0      0;
         0  0 -1  0  0    0  0  1  0  0  1      0;
         0  0  0 -1  0    0  0  0  1  0  0      0;
         0  0  0  0 -1    0  0  0  0  1  0      0];
XU = Polyhedron('H', H_xu);

% safe set
S = Polyhedron('A', [-eye(5); ones(1,nx)], 'b', [zeros(nx,1); 10]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = Dyn(A, F, B, XU);

N = 30;

% Compute attractor defining invariant set
tic
X0 = d.win_always_oneshot_small(S, N, 0.01);
toc
