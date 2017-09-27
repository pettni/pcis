function [ X0 ] = pre(dyn, X, rho)
  % Compute set
  % {x : ∀ p ∃ u ∀ d, x(t+1) + B(rho) ⊆ X}
  % using normal intersection/projection method

  if ~isa(dyn, 'Dyn')
    error('dyn must be an instance of Dyn');
  end

  if nargin < 3
    rho = 0;
  end

  if length(rho) == 1
    rho = rho * ones(dyn.nx,1);
  end

  Xb = X - Polyhedron('A', [eye(dyn.nx); -eye(dyn.nx)], 'b', repmat(rho,2,1));

  X0_A = zeros(0, dyn.nx);
  X0_b = ones(0,1);

  for ip=1:max(1, size(dyn.PV,1))
    A_mat_p = dyn.A;
    F_mat_p = zeros(dyn.nx, 1);
    for jp=1:dyn.np
      A_mat_p = A_mat_p + dyn.Ap{jp} * dyn.PV(ip, jp);
      F_mat_p = F_mat_p + dyn.Fp{jp} * dyn.PV(ip, jp);
    end

    Xd_A = zeros(0, dyn.nx+dyn.nu);
    Xd_b = zeros(0, 1);

    for id=1:max(1, size(dyn.DV,1))
      A_mat_pd = A_mat_p;
      F_mat_pd = F_mat_p;
      for jd=1:dyn.nd
        A_mat_pd = A_mat_pd + dyn.Ad{jd} * dyn.DV(id, jd);
        F_mat_pd = F_mat_pd + dyn.Fd{jd} * dyn.DV(id, jd);
      end
      Xd_A = [Xd_A; Xb.A*[A_mat_pd dyn.B]];
      Xd_b = [Xd_b; Xb.b-Xb.A*F_mat_pd];
    end
    pre_proj = Polyhedron('A', [Xd_A; dyn.XU.A], ...
                          'b', [Xd_b; dyn.XU.b]);
    
    proj = projection(pre_proj, 1:dyn.nx);
    proj.minHRep;

    X0_A = [X0_A; proj.A];
    X0_b = [X0_b; proj.b];
  end
  X0 = Polyhedron('H', [X0_A X0_b]);
  X0.minHRep;
end
