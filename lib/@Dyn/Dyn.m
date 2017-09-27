classdef Dyn
  % Dyn(Ap, Fp, P, B, XU, Ad, Fd, D): 
  %  
  % Create discrete-time system of the form
  %
  %  x(k+1) = (A + ∑d_i Ad{i} + ∑p_i Ap{i}) x(k) + B u + ∑d_i Fd{i} + ∑p_i Fp{i}
  %
  % (x(k),u(k)) ∈ XU input
  % p ∈ P measurable disturbance
  % d ∈ D non-measurable disturbance
  %
  properties (SetAccess=protected)
    A;
    B;
    Ap;
    Fp;
    P;
    XU;
    Ad;
    Fd;
    D;
  end

  methods
    % Constructor
    function d = Dyn(A, B, XU, Ap, Fp, P, Ad, Fd, D)

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

      assert(length(d.Ap) == d.P.Dim)
      assert(length(d.Fp) == d.P.Dim)
      assert(length(d.Ad) == d.D.Dim)
      assert(length(d.Fd) == d.D.Dim)
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

  end
end

