% Graphical solution of the time-optimal control problem of
% reaching the maximal controlled invariant set contained in { (x,v) : x \in [-1, 1] }
% for the double-integrator system

% System matrices, state = [x; dx]
A = [0 1;
     0 0];
B = [0;1];

% Bounds
xmax = 1;
umax = 1;

% Time discretization
dt = 0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Adt = expm(A*dt);
Bdt = dt*B;

X = Polyhedron('A', [1 0; -1 0], 'b', [xmax; xmax])

XU = Polyhedron('H', [0 0 1 umax;
				      0 0 -1 umax]);

d = Dyn(Adt, [], Bdt, XU);

% Invariance
Xinv = win_always(d, X, 0.0001, false, 1);
Xr = [Xinv];

% Reach
for i=1:50
	Xr = [d.pre_proj(Xr(1)) Xr];
end

plot(Xr)
xlim([-2.5 2.5])
ylim([-2.5 2.5])

xlabel('$x$', 'interpreter', 'latex')
ylabel('$\dot x$', 'interpreter', 'latex')
