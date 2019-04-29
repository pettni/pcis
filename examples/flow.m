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
S = Polyhedron('A', [-eye(nx); ones(1,nx)], 'b', [zeros(nx,1); 20]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

d = Dyn(A, F, B, XU);

N = 4;

% Compute attractor defining invariant set
X0 = dyn.win_always_oneshot(S, N, 0.01);

clf; hold on
X = X0;
plot(projection(X0, 3:5), 'alpha', 0.7, 'linestyle', 'none')
for i=1:N
	X = intersect(S, d.pre(X));
	plot(projection(X, 3:5), 'alpha', 0.1)
end
