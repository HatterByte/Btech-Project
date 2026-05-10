function Zr = fb_lin(u, int_u, dn, n, dTe, Te, dTl, Tl, Tl0, Te0, Gr1)
%FB_LIN  Feedback linearization transform — exactly from paper's fb_lin.m
%
%   Confirmed from paper's Simulink source (image):
%     Tl0=314; Te0=290; dac=(-4*dn)*10^(-5); ac=(-4*n-17.3)*10^(-5);
%     Gr=0.0145; Gr1=1.2*Gr;
%     Zr = (1/n^2)*(u*n - dn*int_u) ...
%          - (0.5*(1/Gr1)*(dac*(Tl-Tl0)) + ac*dTl) ...
%          - (0.5*(1/Gr1)*(dac*(Te-Te0)) + ac*dTe);
%
%   INPUTS:
%     u     : control input (rod velocity, δk/k/s)
%     int_u : integral of u (= rho_rod / Gr ≈ cumulative rod movement)
%     dn    : dn_r/dt (rate of change of neutron power)
%     n     : n_r (current neutron power)
%     dTe   : dTe/dt (rate of change of inlet temperature)
%     Te    : current inlet temperature (°C)
%     dTl   : dTl/dt (rate of change of coolant temperature)
%     Tl    : current coolant temperature (°C)
%     Tl0   : reference coolant temperature (314°C)
%     Te0   : reference inlet temperature  (290°C)
%     Gr1   : rod gain with uncertainty (= 1.2*Gr = 1.2*0.0145)

% Rate of change of coolant feedback coefficient
dac = (-4 * dn) * 1e-5;

% Current coolant feedback coefficient
ac = (-4*n - 17.3) * 1e-5;

% Feedback linearization output Zr
% This transforms the nonlinear reactor dynamics into a linear form
% so the ESO can treat it as a simple double integrator.
Zr = (1/n^2) * (u*n - dn*int_u) ...
   - (0.5*(1/Gr1) * (dac*(Tl - Tl0)) + ac*dTl) ...
   - (0.5*(1/Gr1) * (dac*(Te - Te0)) + ac*dTe);

end
