%% Main function to generate tests
function tests = exampleTest
  tests = functiontests(localfunctions);
end

%% Test Functions
function testAlways1(testCase)
  d = Dyn(0.99*[cosd(10) sind(10); -sind(10) cosd(10)]);
  X0 = Polyhedron([eye(2); -eye(2)], ones(4,1));
  W = d.win_always(X0, 0.001);
  cc = W.chebyCenter;
  testCase.assertEqual(cc.r, 1, 'RelTol', 1e-5)
end

function testAlways2(testCase)
  XU = Polyhedron('H', [0 0 1 0 0.1; 0 0 -1 0 0.1; 0 0 0 1 0.1; 0 0 0 -1 0.1]);
  d = Dyn([1 0.2; 0.2 1], [], eye(2), XU);
  X0 = Polyhedron([eye(2); -eye(2)], ones(4,1));
  W = d.win_always(X0, 0.001);
  testCase.assertThat(W <= 1.01*Polyhedron('V', [1 0; 0 1; -1 1; -1 0; 0 -1; 1 -1]), ...
                      matlab.unittest.constraints.IsTrue)
  testCase.assertThat(W >= 0.99*Polyhedron('V', [1 0; 0 1; -1 1; -1 0; 0 -1; 1 -1]), ...
                      matlab.unittest.constraints.IsTrue)

end