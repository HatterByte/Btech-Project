function dx = pwr_nonlinear_dynamics(x, rho, Te)
%PWR_NONLINEAR_DYNAMICS Nonlinear dynamics of the PWR plant.
%
%   INPUTS:
%     x   : State vector [n; c; Tf; Tl]
%     rho : Total reactivity
%     Te  : Coolant inlet temperature (can be a constant or disturbance)
%
%   STATES:
%     n  : Relative neutron power
%     c  : Relative precursor concentration
%     Tf : Average fuel temperature
%     Tl : Average coolant temperature

n  = x(1);
c  = x(2);
Tf = x(3);
Tl = x(4);

lmb = 0.15;
B   = 0.006019;
L   = 0.0001;
ff  = 0.92;

% Variable parameters based on power level
M   = 28*n + 74;
omg = (5/3)*n + 4.933;

uf  = 26.3;
uc  = (160/9)*n + 54.022;

P = 2500;

% Uncertain/modified parameters as per paper
B1   = 0.8*B;
lmb1 = 0.8*lmb;

% Nonlinear ODEs
dn  = (n*(rho - B1)/L) + (B1*c/L);

dc  = lmb1*(n - c);

dTf = (1/uf)*(ff*P*n - omg*Tf + 0.5*omg*(Tl + Te));

dTl = (1/uc)*((1-ff)*P*n + omg*Tf ...
       - 0.5*(2*M + omg)*Tl ...
       + 0.5*(2*M - omg)*Te);

dx = [dn; dc; dTf; dTl];
end
