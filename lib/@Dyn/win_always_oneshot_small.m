function [result] = win_always_oneshot_small(dyn, X, N, eps)
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

  H = [eye(dyn.nx); -eye(dyn.nx)];
  h = [X.chebyCenter.x + eps*ones(dyn.nx,1); eps*ones(dyn.nx,1) - X.chebyCenter.x];

  P = Polyhedron(H, h);

  F = X.A;
  f = X.b;

  Gx = dyn.XU.A(:, 1:dyn.nx);
  Gu = dyn.XU.A(:, dyn.nx+1 : dyn.nx+dyn.nu);
  g = dyn.XU.b;

  Gbar_x_1 = H * dyn.A^N;

  Gbar_x_2_cell = cell([N+1, 1]);
  for i=0:N
    Gbar_x_2_cell{i+1, 1} = F * dyn.A^i;
  end

  Gbar_x_3_cell = cell([N, 1]);  
  for i=0:N-1
    Gbar_x_3_cell{i+1, 1} = Gx * dyn.A^i;
  end

  Gbar_x = [Gbar_x_1; 
            cell2mat(Gbar_x_2_cell); 
            cell2mat(Gbar_x_3_cell)];

  Gbar_u_1_cell = cell([1, N]);
  for i=1:N
    Gbar_u_1_cell{1, i} = H * dyn.A^(N-i) * dyn.B;
  end 

  Gbar_u_2_cell = cell([N+1, N]);
  for j=1:N+1
    for i=1:N
      if j > i
        Gbar_u_2_cell{j, i} = F * dyn.A^(j-i-1) * dyn.B;
      else
        Gbar_u_2_cell{j, i} = zeros(size(F,1), size(dyn.B,2));
      end
    end
  end 

  Gbar_u_3_cell = cell([N, N]);
  for j=1:N
    for i=1:N
      if j==i
        Gbar_u_3_cell{j, i} = Gu;
      elseif j > i
        Gbar_u_3_cell{j, i} = Gx * dyn.A^(j-i-1) * dyn.B;
      else
        Gbar_u_3_cell{j, i} = zeros(size(Gx,1), size(dyn.B,2));
      end
    end
  end 

  K_cell = cell([1, N+1]);
  K_cell{1,1} = zeros(dyn.nx,1);
  for i=2:N+1
    K_cell{1,i} = dyn.A * K_cell{1,i-1} + dyn.F;
  end

  F_K = reshape(F * cell2mat(K_cell), [], 1);
  Gx_K = reshape(Gx * cell2mat(K_cell), [], 1);

  gtilde = [h; zeros(size(g,1) * N, 1); zeros(size(f,1) * (N+1), 1)];

  ghat = [-H*K_cell{1,N}; repmat(f, N+1, 1) - F_K; repmat(g, N, 1) - Gx_K(1:size(Gx,1)*N, :)];

  % gtilde = [h-H*K_cell{1,N}; 
  %           -F_K; 
  %           -Gx_K(1:size(Gx,1)*N, :)];
  % ghat = [zeros(size(h)); repmat(f, N+1, 1); repmat(g, N, 1)];

  Gbar_u = [cell2mat(Gbar_u_1_cell);
            cell2mat(Gbar_u_2_cell);
            cell2mat(Gbar_u_3_cell);];

  sizeT = [size(Gbar_x, 1), size(H, 1)];
  sizeM = [size(Gbar_u, 2), size(H, 2)];

  sizeTvec = sizeT(1) * sizeT(2);
  sizeMvec = sizeM(1) * sizeM(2);

  % T * H = Gbar_x + Gbar_u * M  <=> c1_T * vec(T) = vec(Gbar_x) + c1_M vec(M)
  % for c1_T = H' ⊗ I
  %     c1_M = I ⊗ Gbar_u

  c1_T = kron(H', eye(sizeT(1)));
  c1_M = kron(eye(sizeM(2)), Gbar_u);
  c1_g = zeros(size(c1_T, 1), 1);
  c1_b = vec(Gbar_x);

  % T * h <= gamma * gbar + gtilde  <=> c2_T vec(T) <= gamma * gbar + gtilde
  % for c2_T = h' ⊗ I

  c2_g = -ghat;  
  c2_T = kron(h', eye(sizeT(1)));
  c2_M = zeros(size(c2_T, 1), sizeMvec);
  c2_b = gtilde;

  Aeq = [c1_g c1_T -c1_M];
  beq = [c1_b];

  Aiq = [c2_g c2_T c2_M; 
         -eye(sizeTvec+1) zeros(sizeTvec+1, sizeMvec)];
  biq = [c2_b; zeros(sizeTvec+1, 1)];

  f = [1 zeros(1, sizeTvec + sizeMvec)];

  tic
  [x, fval, exitflag] = linprog(f, Aiq, biq, Aeq, beq);
  toc

  disp("Finishded with exitflag", num2str(exitflag)) 
  exitflag

  result = (1/fval) * P;