%% Main function to generate tests
function tests = exampleTest
  tests = functiontests(localfunctions);
end

%% Test Functions
function testPre(testCase)
  % Identity systems
  A = eye(2);
  dyn = Dyn(A);

  X = Polyhedron('A', [eye(2); -eye(2)], 'b', ones(4,1));
  testCase.assertThat(dyn.pre(X) == X, matlab.unittest.constraints.IsTrue);

  dyn.A = -dyn.A;
  testCase.assertThat(dyn.pre(X) == X, matlab.unittest.constraints.IsTrue);

  % With drift
  dyn = Dyn(A, [1; 0]);
  testCase.assertThat(dyn.pre(X) == Polyhedron('A', [eye(2);-eye(2)], 'b', [0 1 2 1]), matlab.unittest.constraints.IsTrue);

  % With control
  XU = Polyhedron('H', [0 0 1 0 0.1; 0 0 0 1 0.1; 0 0 -1 0 0.1; 0 0 0 -1 0.1]);
  dyn = Dyn(A, [], eye(2), XU);

  testCase.assertThat(dyn.pre(X) == Polyhedron([eye(2); -eye(2)], 1.1*ones(4,1)), ...
                      matlab.unittest.constraints.IsTrue);

  % With control and non-measurable disturbance
  dyn.Ad = {zeros(2,2)};
  dyn.Fd = {[1;0]};
  dyn.D = Polyhedron('V', [0.1; -0.1]);

  testCase.assertThat(dyn.pre(X) == Polyhedron([eye(2); -eye(2)], [1; 1.1; 1; 1.1]), ...
                      matlab.unittest.constraints.IsTrue);

  % With 1D measurable disturbance
  XU = Polyhedron('H', [0 0 1 0 1; 0 0 0 1 1; 0 0 -1 0 1; 0 0 0 -1 1]);
  dyn = Dyn(A, [], eye(2), XU, ...
             {zeros(2,2)}, {[1;0]}, ...
             Polyhedron('V', [1.1; -1.1]));

  testCase.assertThat(dyn.pre(X) == Polyhedron([eye(2); -eye(2)], [0.9; 2; 0.9; 2]), ...
                       matlab.unittest.constraints.IsTrue);

  % With 2D measurable disturbance
  XU = Polyhedron('H', [0 0 1 0 1; 0 0 0 1 1; 0 0 -1 0 1; 0 0 0 -1 1]);
  dyn6 = Dyn(A, [], eye(2), XU, ...
             {zeros(2,2), zeros(2,2)}, {[1;0], [0;1]}, ...
             Polyhedron('V', [0.85 0.85; -0.85 0.85; 0.85 -0.85; -0.85 -0.85]), ...
             [], [], []);
  X0 = dyn6.pre(X);
  testCase.assertThat(X0 == Polyhedron([eye(2); -eye(2)], [1.15; 1.15; 1.15; 1.15]), ...
                       matlab.unittest.constraints.IsTrue);

  % Same with nonmeasurable disturbance
  dyn7 = Dyn(A, [], eye(2), XU, ...
             [], [], [], ...
             {zeros(2,2), zeros(2,2)}, {[1;0], [0;1]}, ...
             Polyhedron('V', [1.1 -1.1; -1.1 1.1; 1.1 1.1; -1.1 -1.1]));

  testCase.assertThat(dyn7.pre(X) == Polyhedron(zeros(0,2), zeros(0,1)), ...
                      matlab.unittest.constraints.IsTrue);

  % Both types
  dyn8 = Dyn(eye(2), [], eye(2), Polyhedron('A', [zeros(2,2) eye(2); zeros(2,2) -eye(2)], 'b', 0.5*ones(4,1)), ...
             {zeros(2,2)}, {[1;0]}, Polyhedron('V', [1.1; -1.1]), ...
             {zeros(2,2)}, {[0;1]}, Polyhedron('V', [0.9; -0.9]));
  testCase.assertThat(dyn8.pre(X) == Polyhedron([eye(2); -eye(2)], [.4; .6; .4; .6]), ...
                      matlab.unittest.constraints.IsTrue);

  % State-dependent non-measurable disturbance
  dyn9 = Dyn(eye(2));
  dyn9.Ew = [0;1];
  dyn9.XW_V = {[0.5 0 0.1], [-0.5 0 -0.1]};
  X = Polyhedron([eye(2); -eye(2)], [4; 2; 0; 2]);
  testCase.assertThat(dyn9.pre(X) == Polyhedron('V', [0 1.9; 0 -1.9; 3.8 0]), ...
                      matlab.unittest.constraints.IsTrue);

  % State-dependent measurable disturbance
  dyn10 = Dyn(eye(2));
  dyn10.Ev = [0;1];
  dyn10.XV_V = {[0.5 0 0.1], [-0.5 0 -0.1]};
  X = Polyhedron([eye(2); -eye(2)], [4; 2; 0; 2]);
  testCase.assertThat(dyn10.pre(X) == Polyhedron('V', [0 1.9; 0 -1.9; 3.8 0]), ...
                      matlab.unittest.constraints.IsTrue);

end

%% Test Functions
function testPreRc(testCase)
  % Identity system
  A = eye(2);
  dyn1 = Dyn(A);
  dyn1p = Dyn(-A);

  X = Polyhedron('A', [eye(2); -eye(2)], 'b', ones(4,1));

  testCase.assertThat(dyn1.pre_rc(X) == X, matlab.unittest.constraints.IsTrue);
  testCase.assertThat(dyn1p.pre_rc(X) == X, matlab.unittest.constraints.IsTrue);

  % With drift
  dyn = Dyn(A, [1; 0]);
  testCase.assertThat(dyn.pre(X) == Polyhedron('A', [eye(2);-eye(2)], 'b', [0 1 2 1]), matlab.unittest.constraints.IsTrue);

  % With control
  XU = Polyhedron('H', [0 0 1 0 0.1; 0 0 0 1 0.1; 0 0 -1 0 0.1; 0 0 0 -1 0.1]);
  dyn3 = Dyn(A, [], eye(2), XU);

  testCase.assertThat(dyn3.pre_rc(X) == Polyhedron([eye(2); -eye(2)], 1.1*ones(4,1)), ...
                      matlab.unittest.constraints.IsTrue);

  % With control and non-measurable disturbance
  dyn4 = Dyn(A, [], eye(2), XU, [], [], [],...
             {zeros(2,2)}, {[1;0]}, Polyhedron('H', [1 0.1; -1 0.1]));

  testCase.assertThat(dyn4.pre_rc(X) == Polyhedron([eye(2); -eye(2)], [1; 1.1; 1; 1.1]), ...
                      matlab.unittest.constraints.IsTrue);

  % With 1D measurable disturbance
  XU = Polyhedron('H', [0 0 1 0 1; 0 0 0 1 1; 0 0 -1 0 1; 0 0 0 -1 1]);
  dyn5 = Dyn(A, [], eye(2), XU, ...
             {zeros(2,2)}, {[1;0]}, ...
             Polyhedron('H', [1 1.1; -1 1.1]), ...
             [], [], []);

  testCase.assertThat(dyn5.pre_rc(X) == Polyhedron([eye(2); -eye(2)], [0.9; 2; 0.9; 2]), ...
                       matlab.unittest.constraints.IsTrue);

  % With 2D measurable disturbance
  XU = Polyhedron('H', [0 0 1 0 1; 0 0 0 1 1; 0 0 -1 0 1; 0 0 0 -1 1]);
  dyn6 = Dyn(A, [], eye(2), XU, ...
             {zeros(2,2), zeros(2,2)}, {[1;0], [0;1]}, ...
             Polyhedron('A', [eye(2); -eye(2)], 'b', 0.85*ones(4,1)), ...
             [], [], []);
  X0 = dyn6.pre_rc(X);
  testCase.assertThat(X0 == Polyhedron([eye(2); -eye(2)], [1.15; 1.15; 1.15; 1.15]), ...
                       matlab.unittest.constraints.IsTrue);

  % Same with nonmeasurable disturbance
  dyn7 = Dyn(A, [], eye(2), XU, ...
             [], [], [], ...
             {zeros(2,2), zeros(2,2)}, {[1;0], [0;1]}, ...
             Polyhedron('A', [eye(2); -eye(2)], 'b', 1.1*ones(4,1)));

  testCase.assertThat(dyn7.pre_rc(X) == Polyhedron(zeros(0,2), zeros(0,1)), ...
                      matlab.unittest.constraints.IsTrue);

  % Both types
  dyn8 = Dyn(eye(2), [], eye(2), Polyhedron('A', [zeros(2,2) eye(2); zeros(2,2) -eye(2)], 'b', 0.5*ones(4,1)), ...
             {zeros(2,2)}, {[1;0]}, Polyhedron('H', [1 1.1; -1 1.1]), ...
             {zeros(2,2)}, {[0;1]}, Polyhedron('H', [1 0.9; -1 0.9]));
  testCase.assertThat(dyn8.pre_rc(X) == Polyhedron([eye(2); -eye(2)], [.4; .6; .4; .6]), ...
                      matlab.unittest.constraints.IsTrue);

end
