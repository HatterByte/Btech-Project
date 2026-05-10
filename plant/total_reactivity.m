function rho = total_reactivity(n, Tf, Tl, Te, rho_rod)
%TOTAL_REACTIVITY Computes the total reactivity including feedback.
%
%   INPUTS:
%     n       : Relative neutron power
%     Tf      : Fuel temperature
%     Tl      : Coolant temperature
%     Te      : Coolant inlet temperature
%     rho_rod : Reactivity from control rods

Tf0 = 650;
Tl0 = 314;
Te0 = 290;

% Temperature feedback coefficients
af = (n - 4.24)*1e-5;
ac = (-4*n - 17.3)*1e-5;

% Total reactivity formula
rho = rho_rod ...
    + af*(Tf - Tf0) ...
    + 0.5*ac*(Tl - Tl0) ...
    + 0.5*ac*(Te - Te0);

end
