% Finding a controlled invariant set via a one-shot method

% System matrices, state = [x; dx; ddx]
A = [0 1 0;
     0 0 1;
     0 0 0];
B = [0;0;1];
F = [0.1;0;0];

% Bounds
xmax = 1;
dxmax = 1;
ddxmax = 1;
umax = 1;

% Time discretization
dt = 0.25;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Adt = expm(A*dt);
Bdt = dt*B;
Fdt = dt*F;

S = Polyhedron('A', [eye(3); -eye(3)], 'b', [xmax; dxmax; ddxmax; xmax; dxmax; ddxmax]);

XU = Polyhedron('H', [0 0 0 1 umax;
				      0 0 0 -1 umax]);

d = Dyn(Adt, Fdt, Bdt, XU);

N = 20;

% Compute attractor defining invariant set
% tic
X0 = d.win_always_oneshot(S, N, 0.1);
% toc

% Now X0 is contained in pre^N(X0, S) for pre(X, S) := S ∩ pre(X)

% Plot the iterations
clf; hold on
X = X0;
for i=1:N
	X = intersect(S, d.pre(X));
	plot(X, 'alpha', 0.1, 'linestyle', 'none')
end

% Plot auxiliary sets
plot(X0, 'color', 'blue', 'alpha', 0.5)
plot(S, 'edgecolor', 'blue', 'alpha', 0)

axis equal