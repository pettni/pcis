%% poly_dist: compute radius of largest ball around x contained in P
function [r, grad] = poly_ball(P, x)

	A = P.A;
	b = P.b;

    n = sqrt(sum(A.*A,2));

	% normalize 0'*x<=+/-b to 0'*x<= +/- sign(b)
	% (correct sign is important as not to mess with trivial infeasibility)
	ZeroRows = (n<1e-10);
	n(ZeroRows) = 1;
	b(ZeroRows) = sign(b(ZeroRows));

	% normalize each half-space (0*x<=b will be left intact)
	nA = A ./ repmat(n,1,size(A,2));
	nb = b ./ n;

	Aiq = ones(size(nA, 1), 1);
	biq = nb - nA * x;

	[x, fval, exitflag,output,lambda] = linprog([-1], Aiq, biq);

	r = -fval;
	grad = -nA(find(lambda.ineqlin, 1, 'first'), :);