% Finding a controlled invariant set via a one-shot method

% System matrices, state = [x; dx]
A = [0 1;
     0 0];
B = [0;1];
F = [0.1;0];
% Bounds
xmax = 1;
umax = 1;

% Time discretization
dt = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Adt = expm(A*dt);
Bdt = dt*B;
Fdt = dt*F;

XU = Polyhedron('H', [0 0 1 umax;
				      0 0 -1 umax]);

d = Dyn(Adt, Fdt, Bdt, XU);

S = Polyhedron('A', [eye(2); -eye(2)], 'b', [xmax; xmax; xmax; xmax]);

N = 20;

% Compute attractor defining invariant set
X0 = d.win_always_oneshot(S, N, 0.15);

% Now X0 should be contained in pre^N(X0, S) for pre(X, S) := S ∩ pre(X)

% Plot the iterations
clf; hold on
X = X0;
for i=1:N
	X = intersect(S, d.pre(X));
	plot(X, 'alpha', 0.1)
end

% Plot auxiliary sets
plot(X0, 'linestyle', '--', 'alpha', 0)
plot(S, 'edgecolor', 'blue', 'alpha', 0)

% Compare with maximal
plot(d.win_always(S), 'edgecolor', 'green', 'alpha', 0)

axis equal