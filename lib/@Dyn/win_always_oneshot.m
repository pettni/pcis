function [result] = win_always_oneshot(dyn, X, N, eps)
% One-shot search for maximal invariant set contained in X 
% 
% by looking for a set P such that P
% P ⊆ pre^N(P, X),    where pre(P, X) = X ∩ pre(P)
% which is sufficient for existence of a control-invariant set
%
% Reference: Fiacchini, M., & Alamir, M. (2018). Computing control 
%            invariant sets in high dimension is easy.
%            http://arxiv.org/abs/1810.10372

  if nargin < 4
  	eps = 1e-3;
  end

  P = Polyhedron('A', [eye(dyn.nx); -eye(dyn.nx)], 'b', eps*ones(2*dyn.nx, 1));

  if dyn.is_statedep_input()
  	error("state-dependent input bounds not implemented");
  end

  if dyn.np > 0 || dyn.nv > 0
  	error("method does not support parameters");
  end

  if dyn.nd > 0 || dyn.nw > 0
  	error("method does not support disturbance");
  end

  if norm(dyn.F) > 0
  	error("drift not implemented");
  end

  if not(X.contains(zeros(dyn.nx,1)))
  	error("0 must be contained in X")
  end

  H = P.A;
  h = P.b;

  G = dyn.XU.A(:, dyn.nx+1 : dyn.nx+dyn.nu);
  g = dyn.XU.b;

  Gbar_diag = repmat({G}, 1, N);
  Gbar = blkdiag(H * dyn.A^N, Gbar_diag{:});
  for i=1:N
  	Gbar(1:size(H, 1), dyn.nx + dyn.nu*(i-1) + 1: dyn.nx + dyn.nu*i) ...
  						= H * dyn.A^(i-1) * dyn.B;
  end

  Hbar = [H zeros(size(H,1), dyn.nu * N)];

  gtilde = [h; zeros(size(g,1) * N, 1)];
  ghat = [zeros(size(h,1), 1); repmat(g, N, 1)];

  if nargin == 4
  	% add intersection condition
  	F = X.A;
  	f = X.b;

  	Gbar_l_cell = cell([N+1, N+1]);

  	% left column
  	for i=1:N+1
  		Gbar_l_cell{i, 1} = F * dyn.A^(N-i+1);
  	end

  	for j=2:N+1  % for each column
  		for i=1:N+1   % for each row
  			if i < j
	  			Gbar_l_cell{i, j} = F * dyn.A^(j-i-1) * dyn.B;
	  		else
	  			Gbar_l_cell{i, j} = zeros(size(F,1), dyn.nu);
	  		end
  		end
  	end

	Gbar = [Gbar; cell2mat(Gbar_l_cell)];

	ghat = [ghat; repmat(f, N+1, 1)];
	gtilde = [gtilde; zeros(size(f,1) * (N+1), 1)];
  end

  Q = [eye(dyn.nx) repmat(zeros(dyn.nx, dyn.nu), 1, N)];

  sizeT = [size(Gbar, 1), size(Hbar, 1)];
  sizeM = [size(Gbar, 2), size(Hbar, 2)];

  sizeTvec = sizeT(1) * sizeT(2);
  sizeMvec = sizeM(1) * sizeM(2);

  % T * Hbar = Gbar * M  <=> c1_T * vec(T) = c1_M vec(M)
  % for c1_T = Hbar' ⊗ I
  %     c1_M = I ⊗ Gbar

  c1_T = kron(Hbar', eye(sizeT(1)));
  c1_M = kron(eye(sizeM(2)), Gbar);
  c1_g = zeros(size(c1_T, 1), 1);
  c1_b = zeros(size(c1_T, 1), 1);

  % T * h <= gbar  <=> c2_T vec(T) <= gamma * gbar + gtilde
  % for c2_T = h' ⊗ I

  c2_g = -ghat;  
  c2_T = kron(h', eye(sizeT(1)));
  c2_M = zeros(size(c2_T, 1), sizeMvec);
  c2_b = gtilde;

  % Q = Q M  <=> vec(Q) = c3_M * vec(M)
  % for c3_M = I ⊗ Q

  c3_M = kron(eye(sizeM(2)), Q);
  c3_T = zeros(size(c3_M, 1), sizeTvec);
  c3_g = zeros(size(c3_M, 1), 1);
  c3_b = vec(Q);

  Aeq = [c1_g c1_T -c1_M; 
         c3_g c3_T -c3_M];
  beq = [c1_b; c3_b];

  Aiq = [c2_g c2_T c2_M; 
         -eye(sizeTvec+1) zeros(sizeTvec+1, sizeMvec)];
  biq = [c2_b; zeros(sizeTvec+1, 1)];

  f = [1 zeros(1, sizeTvec + sizeMvec)];

  [x, fval, exitflag] = linprog(f, Aiq, biq, Aeq, beq);

  result = (1/fval) * P;