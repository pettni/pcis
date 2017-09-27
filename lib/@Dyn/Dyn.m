classdef Dyn
  % Dyn(A, B, XU, Ap, Fp, P, Ad, Fd, D, E, XW_V): 
  % Discrete-time system of the form
  %
  %  x(k+1) = (A + ∑d_i Ad{i} + ∑p_i Ap{i}) x(k) + B u(k) ...
  %           + ∑d_i(k) Fd{i} + ∑p_i(k) Fp{i} + E w(k)
  %
  % u input such that (x(k),u(k)) ∈ XU
  % p(k) ∈ P measurable disturbance
  % d(k) ∈ D non-measurable disturbance
  % w(k) ∈ conv_i(XW_V{i} [x(k); 1]) non-measurable state-dependent disturbance
  %
  % INPUTS:
  %
  % A: (nx x nx) matrix
  % B: (nx x nu) matrix
  % E: (nx x nw) matrix
  %
  % Ap: cell of np matrices of size (nx x nx)
  % Fp: cell of np matrices of size (nx x 1)  
  % Ad: cell of nd matrices of size (nx x nx)
  % Fd: cell of nd matrices of size (nx x 1) 
  %
  % XU: (nx+nu)-dim Polyhedron
  % P:  np-dim Polyhedron
  % D:  nd-dim Pokyhedron
  %
  % XW_V: cell of matrices of size (nw x nx+1)  
  %
  properties (SetAccess=protected)
    A;
    B;
    XU;
    Ap;
    Fp;
    P;
    Ad;
    Fd;
    D;
    E;
    XW_V;
  end

  methods
    % Constructor
    function d = Dyn(A, B, XU, Ap, Fp, P, Ad, Fd, D, E, XW_V)

      nx = size(A,2);

      d.A = A;

      if nargin < 2 || isempty(B)
        d.B = zeros(nx,0);
        d.XU = Polyhedron('H', [zeros(1,nx) 1]);
      else
        d.B = B;
        d.XU = XU;
      end

      if nargin < 4 || isempty(Ap)
        d.Ap = {};
        d.Fp = {};
        d.P = Polyhedron; 
      else
        d.Ap = Ap;
        d.Fp = Fp;
        d.P = P;
      end

      if nargin < 7 || isempty(Ad)
        d.Ad = {};
        d.Fd = {};
        d.D = Polyhedron; 
      else
        d.Ad = Ad;
        d.Fd = Fd;
        d.D = D;
      end

      if nargin < 10 || isempty(E)
        d.E = zeros(nx,0);
        d.XW_V = {};
      else
        d.E = E;
        d.XW_V = XW_V;
      end

      % Checks
      assert(size(d.A, 1) == size(d.A,2))
      assert(size(d.B, 1) == size(d.A,1))
      assert(size(d.E, 1) == size(d.A,1))

      assert(d.XU.Dim == size(d.A,2) + size(d.B,2))

      assert(length(d.Ap) == d.P.Dim)
      assert(length(d.Fp) == d.P.Dim)
      assert(length(d.Ad) == d.D.Dim)
      assert(length(d.Fd) == d.D.Dim)

      for i=1:d.P.Dim
        assert(size(d.Ap{i}, 1) == size(d.A,1))
        assert(size(d.Ap{i}, 2) == size(d.A,1))
        assert(size(d.Fp{i}, 1) == size(d.A,1))
        assert(size(d.Fp{i}, 2) == 1)
      end

      for i=1:d.D.Dim
        assert(size(d.Ad{i}, 1) == size(d.A,1))
        assert(size(d.Ad{i}, 2) == size(d.A,1))
        assert(size(d.Fd{i}, 1) == size(d.A,1))
        assert(size(d.Fd{i}, 2) == 1)
      end

      for i=1:length(d.XW_V)
        assert(size(d.E,2) == size(d.XW_V{i}, 1))
        assert(size(d.A,2)+1 == size(d.XW_V{i}, 2))
      end

    end

    function n = nx(d)
      n = size(d.A,2);
    end
    function n = nu(d)
      n = size(d.B,2);
    end
    function n = nd(d)
      n = length(d.Ad);
    end
    function n = np(d)
      n = length(d.Ap);
    end
    function n = nw(d)
      n = size(d.E,2);
    end
    function X0 = pre(d, X, rho)
      if nargin < 3
        rho = 0;
      end
      X0 = pre_proj(d, X, rho);
    end
  end
end

