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

  if dyn.np > 0 || dyn.nv > 0
    error("method does not support parameters");
  end

  if dyn.nd > 0 || dyn.nw > 0
    error("method does not support disturbance");
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  mas_A = [X.A zeros(size(X.A,1), dyn.nu);
           dyn.XU.A];
  mas_b = [X.b; 
           dyn.XU.b];

  mas_poly = Polyhedron('A', mas_A, 'b', mas_b, 'Ae', [dyn.A-eye(dyn.nx) dyn.B], 'be', [-dyn.F]);
  mas_poly.chebyCenter.x;
  xc = mas_poly.chebyCenter.x(1:dyn.nx);

  % place equilibrium at origin
  A = dyn.A;
  B = dyn.B;
  K = dyn.F + A * xc - xc;

  % ball around origin is target polyhedron Omega
  H = [eye(dyn.nx); -eye(dyn.nx)];
  h = [eps*ones(2*dyn.nx,1)];
  Omega = Polyhedron(H, h);

  % translate constraints
  X = X - xc;
  F = X.A;
  f = X.b;

  % translate input constraints
  XU = dyn.XU - [xc; zeros(dyn.nu, 1)];
  Gx = XU.A(:, 1:dyn.nx);
  Gu = XU.A(:, dyn.nx+1 : dyn.nx+dyn.nu);
  g = XU.b;

  % check that origin is equilibrium
  uc = mas_poly.chebyCenter.x(dyn.nx+1:end);

  if ~X.contains(Omega)
    error("containment error")
  end

  if norm(B * uc + K) > 1e-6
    error('zero is not equilibrium, problem with translation');
  end

  if ~X.contains(zeros(dyn.nx, 1))
    error('translated X does not contain origin, problem with translation')
  end

  if ~XU.contains([zeros(dyn.nx, 1); uc])
    error('translated XU does not contain equilibrium, problem with translation')
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Gbar_x_2_cell = cell([N+1, 1]);
  for i=0:N
    Gbar_x_2_cell{i+1, 1} = F * A^i;
  end

  Gbar_x_3_cell = cell([N, 1]);  
  for i=0:N-1
    Gbar_x_3_cell{i+1, 1} = Gx * A^i;
  end

  Gbar_x = [H * A^N; 
            cell2mat(Gbar_x_2_cell); 
            cell2mat(Gbar_x_3_cell)];

  Gbar_u_1_cell = cell([1, N]);
  for i=1:N
    Gbar_u_1_cell{1, i} = H * A^(N-i) * B;
  end 

  Gbar_u_2_cell = cell([N+1, N]);
  for j=1:N+1
    for i=1:N
      if j > i
        Gbar_u_2_cell{j, i} = F * A^(j-i-1) * B;
      else
        Gbar_u_2_cell{j, i} = zeros(size(F,1), size(B,2));
      end
    end
  end 

  Gbar_u_3_cell = cell([N, N]);
  for j=1:N
    for i=1:N
      if j==i
        Gbar_u_3_cell{j, i} = Gu;
      elseif j > i
        Gbar_u_3_cell{j, i} = Gx * A^(j-i-1) * B;
      else
        Gbar_u_3_cell{j, i} = zeros(size(Gx,1), size(B,2));
      end
    end
  end 

  K_cell = cell([1, N+1]);
  K_cell{1,1} = zeros(dyn.nx,1);
  for i=2:N+1
    K_cell{1,i} = A * K_cell{1,i-1} + K;
  end

  F_K = reshape(F * cell2mat(K_cell), [], 1);
  Gx_K = reshape(Gx * cell2mat(K_cell), [], 1);

  gtilde = [h; 
            zeros(size(g,1) * N, 1); 
            zeros(size(f,1) * (N+1), 1)];

  ghat = [-H*K_cell{1,N+1};
          repmat(f, N+1, 1) - F_K; 
          repmat(g, N, 1) - Gx_K(1:size(Gx,1)*N, :)];

  Gbar_u = [cell2mat(Gbar_u_1_cell);
            cell2mat(Gbar_u_2_cell);
            cell2mat(Gbar_u_3_cell);];

  sizeT = [size(Gbar_x, 1), size(H, 1)];
  sizeR = [size(Gbar_u, 2), size(H, 2)];

  sizeTvec = sizeT(1) * sizeT(2);
  sizeRvec = sizeR(1) * sizeR(2);

  sizer = N * dyn.nu;

  % T * H = Gbar_x + Gbar_u * R  
  %     <=> c1_T * vec(T) = vec(Gbar_x) + c1_R vec(R)
  % for c1_T = H' ⊗ I
  %     c1_R = I ⊗ Gbar_u

  c1_T = kron(H', eye(sizeT(1)));
  c1_R = -kron(eye(sizeR(2)), Gbar_u);
  c1_g = zeros(size(c1_T, 1), 1);
  c1_r = zeros(size(c1_T, 1), sizer);
  c1_b = vec(Gbar_x);

  % T * h <= gamma * ghat + gtilde - Gbar_u r  
  %     <=>  c2_T vec(T) <= gamma * ghat + gtilde - Gbar_u r
  % for c2_T = h' ⊗ I

  c2_g = -ghat;  
  c2_T = kron(h', eye(sizeT(1)));
  c2_R = zeros(size(c2_T, 1), sizeRvec);
  c2_r = Gbar_u;
  c2_b = gtilde;


  Aeq = [c1_g c1_T c1_R c1_r];
  beq = [c1_b];

  Aiq = [c2_g c2_T c2_R c2_r; 
         -eye(sizeTvec+1) zeros(sizeTvec+1, sizeRvec + sizer)];
  biq = [c2_b; zeros(sizeTvec+1, 1)];

  f = [1 zeros(1, sizeTvec + sizeRvec + sizer)];

  [x, fval, exitflag] = linprog(f, Aiq, biq, Aeq, beq);

  if exitflag < 0
    disp(['Finished with exitflag ', num2str(exitflag), ', probably infeasible']) 
    result = [];
  else
    result = xc + (1/fval) * Omega;
  end