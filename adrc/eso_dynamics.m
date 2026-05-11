function dz = eso_dynamics(z, y, u, varargin)
%ESO_DYNAMICS Linear ESO dynamics.
%
%   Paper-matched form:
%       dz = Af*z + Bf*u + Lf*(y - Cf*z)
%
%   Legacy form:
%       dz = eso_dynamics(z, y, u, wo, b0)

if nargin == 4
    p = varargin{1};
    dz = p.Af*z + p.Bf*u + p.Lf*(y - p.Cf*z);
    return;
end

wo = varargin{1};
b0 = varargin{2};

L1 = 3 * wo;
L2 = 3 * wo^2;
L3 = wo^3;

dz = [z(2) + L1 * (y - z(1));
      z(3) + b0 * u + L2 * (y - z(1));
      L3 * (y - z(1))];
end
