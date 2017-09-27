%% Main function to generate tests
function tests = exampleTest
  tests = functiontests(localfunctions);
end

%% Test Functions
function testPre(testCase)
  % Identity system
  A = eye(2);
  dyn1 = Dyn({A}, {zeros(2,1)});

  X = Polyhedron('A', [eye(2); -eye(2)], 'b', ones(4,1));

  Xp = dyn1.pre(X);
  testCase.assertThat(Xp == X, matlab.unittest.constraints.IsTrue);

  % With drift
  dyn2 = Dyn({A}, {[1; 1]});
  Xp = dyn2.pre(X);
  real_pre = Polyhedron([eye(2); -eye(2)], [0;0;2;2]);

  testCase.assertThat(Xp == real_pre, matlab.unittest.constraints.IsTrue);

  % With control
  XU = Polyhedron('H', [0 0 1 0 0.1; 0 0 0 1 0.1; 0 0 -1 0 0.1; 0 0 0 -1 0.1]);
  dyn3 = Dyn({A}, {[0;0]}, 1, eye(2), XU);
  Xp = dyn3.pre(X);
  real_pre = Polyhedron([eye(2); -eye(2)], 1.1*ones(4,1));

  testCase.assertThat(Xp == real_pre, matlab.unittest.constraints.IsTrue);

  % With control and non-measurable disturbance
  dyn4 = Dyn({A}, {[0;0]}, 1, eye(2), XU, ...
             {zeros(2,2)}, {[1;0]}, [0.1; -0.1]);
  Xp = dyn4.pre(X);
  real_pre = Polyhedron([eye(2); -eye(2)], [1; 1.1; 1; 1.1]);

  testCase.assertThat(Xp == real_pre, matlab.unittest.constraints.IsTrue);

  % With measurable disturbance
  XU = Polyhedron('H', [0 0 1 0 1; 0 0 0 1 1; 0 0 -1 0 1; 0 0 0 -1 1]);
  dyn5 = Dyn({A, zeros(2,2), zeros(2,2)}, {[0;0], [0;1], [1;0]}, [1 1.1 1.1; 1 -1.1 1.1; 1 1.1 -1.1; 1 -1.1 -1.1], eye(2), XU);
  Xp = dyn5.pre(X);
  real_pre = Polyhedron([eye(2); -eye(2)], [0.9; 0.9; 0.9; 0.9]);

  testCase.assertThat(Xp == real_pre, matlab.unittest.constraints.IsTrue);

  % Same with nonmeasurable disturbance
  dyn6 = Dyn({A}, {[0;0]}, [1], eye(2), XU, {zeros(2,2), zeros(2,2)}, ...
             {[1;0], [0;1]}, [1.1 -1.1; -1.1 1.1; 1.1 1.1; -1.1 -1.1]);
  Xp = dyn6.pre(X);
  real_pre = Polyhedron(zeros(0,2), zeros(0,1));

  testCase.assertThat(Xp == real_pre, matlab.unittest.constraints.IsTrue);

end

