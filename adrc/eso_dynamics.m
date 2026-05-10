function dz = eso_dynamics(z, y, u, wo, b0)
%ESO_DYNAMICS Dynamics of the Linear Extended State Observer (LESO).
%
%   Reduced-order model:
%     x1_dot = x2
%     x2_dot = x3 + b0*u
%     x3_dot = h  (unknown disturbance rate)
%
%   Observer Equations:
%     dz1 = z2 + L1*(y - z1)
%     dz2 = z3 + b0*u + L2*(y - z1)
%     dz3 = L3*(y - z1)
%
%   INPUTS:
%     z  : Observer states [z1; z2; z3]
%     y  : Plant output (measured power deviation)
%     u  : Control input
%     wo : Observer bandwidth
%     b0 : Gain parameter

% Gains parameterized by bandwidth wo
L1 = 3 * wo;
L2 = 3 * wo^2;
L3 = wo^3;

z1 = z(1);
z2 = z(2);
z3 = z(3);

dz1 = z2 + L1 * (y - z1);
dz2 = z3 + b0 * u + L2 * (y - z1);
dz3 = L3 * (y - z1);

dz = [dz1; dz2; dz3];
end
