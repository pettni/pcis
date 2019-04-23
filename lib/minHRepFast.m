%% minHrepFast: wrapper for fast minHrep for MPT polyhedron 
function [p] = minHrepFast(p)
	[Amin, bmin] = indicate_nonredundant_halfplanes(p.A, p.b);
	p = Polyhedron(Amin, bmin);
end