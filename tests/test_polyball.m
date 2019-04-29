%% Main function to generate tests
function tests = exampleTest
  tests = functiontests(localfunctions);
end

%% Test Functions
function testPolyBall(testCase)

	P = Polyhedron(eye(2), [1;1]);   % gradient ambiguous (non-smooth)
	[r, grad] = poly_ball(P, [1;1]);
	testCase.assertEqual(r, 0, 'RelTol', 1e-9) 

	[r, grad] = poly_ball(P, [0.2;0.3]);
	testCase.assertEqual(r, 0.7, 'RelTol', 1e-9)
	testCase.assertEqual(grad, [0 -1], 'RelTol', 1e-9)

	[r, grad] = poly_ball(P, [0;0.9]);
	testCase.assertEqual(r, 0.1, 'RelTol', 1e-9)
	testCase.assertEqual(grad, [0 -1], 'RelTol', 1e-9)
end