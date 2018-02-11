% Compute an invariant set for a triple integrator system

% System matrices, state = [x; dx; ddx]
A = [0 1 0;
     0 0 1;
     0 0 0];
B = [0;0;1];

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

X = Polyhedron('A', [eye(3); -eye(3)], 'b', [xmax; dxmax; ddxmax; xmax; dxmax; ddxmax])

XU = Polyhedron('H', [0 0 0 1 umax;
				      0 0 0 -1 umax]);

d = Dyn(Adt, [], Bdt, XU);

Xinv = win_always(d, X, 0.05, false, 1);

plot(Xinv, 'alpha', 0.75)
xlabel('$x$', 'interpreter', 'latex')
ylabel('$\dot x$', 'interpreter', 'latex')
zlabel('$\ddot x$', 'interpreter', 'latex')