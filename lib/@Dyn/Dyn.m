classdef Dyn
  % DYN: Create a Dyn object.
  % ========================================
  %
  % SYNTAX
  % ------
  %
  %   dyn = Dyn(Ap, Fp, P, B, XU, Ad, Fd, D)
  % 
  % DESCRIPTION
  % ------    
  %   Discrete-time system of the form 
  %
  %  x(t+1) = (∑d_i Ad{i} + ∑p_i Ap{i}) x(t) + B u + ∑d_i Fd{i} + ∑p_i Fp{i}
  %
  % (x,u) ∈ XU input
  % p ∈ conv(PV) measurable disturbance
  % d ∈ conv(DV) non-measurable disturbance
  %
  properties (SetAccess=protected)
    Ap;
    Fp;
    PV;
    B;
    XU;
    Ad;
    Fd;
    DV;
  end

  methods
    % Constructor
    function d = Dyn(Ap, Fp, PV, B, XU, Ad, Fd, DV)

      nx = size(Ap{1},2);

      d.Ap = Ap;
      d.Fp = Fp;

      if nargin < 3
        assert length(Ap) == 1
        d.PV = 1 
      else
        d.PV = PV;
      end

      if nargin < 4
        d.B = zeros(nx, 0);
        d.XU = Polyhedron('H', [zeros(1,nx) 1]);
      else
        d.B = B;
        d.XU = XU;
      end

      if nargin < 7
        d.Ad = {zeros(nx, nx)};
        d.Fd = {zeros(nx, 1)};
        d.DV = 1; 
      else
        d.Ad = Ad;
        d.Fd = Fd;
        d.DV = DV;
      end

      assert(length(d.Ap) == size(d.PV,1))
      assert(length(d.Fp) == size(d.PV,1))
      assert(length(d.Ad) == size(d.DV,1))
      assert(length(d.Fd) == size(d.DV,1))
    end

    function n = nx(d)
      n = size(d.Ap{1},2);
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

