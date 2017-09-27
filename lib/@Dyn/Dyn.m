classdef Dyn
  % Dyn(Ap, Fp, P, B, XU, Ad, Fd, D): 
  %  
  % Create discrete-time system of the form
  %
  %  x(k+1) = (∑d_i Ad{i} + ∑p_i Ap{i}) x(k) + B u + ∑d_i Fd{i} + ∑p_i Fp{i}
  %
  % (x(k),u(k)) ∈ XU input
  % p ∈ conv(PV) measurable disturbance
  % d ∈ conv(DV) non-measurable disturbance
  %
  properties (SetAccess=protected)
    A;
    B;
    Ap;
    Fp;
    PV;
    XU;
    Ad;
    Fd;
    DV;
  end

  methods
    % Constructor
    function d = Dyn(A, B, XU, Ap, Fp, PV, Ad, Fd, DV)

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
        d.PV = []; 
      else
        d.Ap = Ap;
        d.Fp = Fp;
        d.PV = PV;
      end

      if nargin < 7 || isempty(Ad)
        d.Ad = {};
        d.Fd = {};
        d.DV = []; 
      else
        d.Ad = Ad;
        d.Fd = Fd;
        d.DV = DV;
      end

      assert(length(d.Ap) == size(d.PV,2))
      assert(length(d.Fp) == size(d.PV,2))
      assert(length(d.Ad) == size(d.DV,2))
      assert(length(d.Fd) == size(d.DV,2))
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

