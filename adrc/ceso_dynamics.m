function [dz, x_hat] = ceso_dynamics(z, y, u, p)
%CESO_DYNAMICS Two-level cascade ESO used in CESO-ADRC.
%
%   z = [z_level_1; z_level_2], each level is 3x1.

z1 = z(1:3);
z2 = z(4:6);

y0 = y;
y1 = p.Cf * z1;

dz1 = p.Af*z1 + p.Bf*u + p.L1*(y0 - p.Cf*z1);
dz2 = p.Af*z2 + p.Bf*u + p.L2*(y1 - p.Cf*z2) + p.Gamma*z1;

% CESO estimate used by the ADRC control law.
x_hat = [z2(1); z2(2); z2(3) + z1(3)];

dz = [dz1; dz2];
end
