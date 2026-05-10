function dx = pwr_plant_paper(x, rho_rod, Te, d_rho, Tf_ref, Tl_ref, Te_ref)
%PWR_PLANT_PAPER  Nonlinear PWR plant dynamics.
%
%   Exactly matches npp_nl.m and tot_react.m from paper's Simulink source.
%
%   Equations confirmed from paper images:
%     tot_react: rho = rho_rod + af*(Tf-Tf0) + 0.5*ac*(Tl-Tl0) + 0.5*ac*(Te-Te0)
%                af = (n-4.24)*1e-5,  ac = (-4*n-17.3)*1e-5
%                Tf0=650, Tl0=314, Te0=290
%
%     npp_nl:   dn  = (n*(rho-B1)/l) + (B1*c/l)
%               dc  = lmb1*(n-c)
%               dTf = (1/uf)*(ff*P*n - omg*Tf + 0.5*omg*(Tl+Te))
%               dTl = (1/uc)*((1-ff)*P*n + omg*Tf - 0.5*(2M+omg)*Tl + 0.5*(2M-omg)*Te)
%
%   States x (4x1): [n; c; Tf; Tl]

% Constants (from npp_nl.m)
lmb = 0.15;
B   = 0.006019;
l   = 1e-4;
ff  = 0.92;
P   = 2500;
uf  = 26.3;

% Uncertain parameters (paper uses 0.8× nominal in controller model)
lmb1 = 0.8 * lmb;
B1   = 0.8 * B;

% States
n  = x(1);
c  = x(2);
Tf = x(3);
Tl = x(4);

% State-dependent parameters (nonlinear)
M   = 28*n + 74;
omg = (5/3)*n + 4.933;
uc  = (160/9)*n + 54.022;

% Feedback coefficients (from tot_react.m, confirmed)
af = (n - 4.24) * 1e-5;
ac = (-4*n - 17.3) * 1e-5;

% Total reactivity (from tot_react.m)
rho = rho_rod ...
    + af*(Tf - Tf_ref) ...
    + 0.5*ac*(Tl - Tl_ref) ...
    + 0.5*ac*(Te - Te_ref) ...
    + d_rho;

% Point kinetics (from npp_nl.m)
dn  = (n*(rho - B1)/l) + (B1*c/l);
dc  = lmb1*(n - c);

% Heat exchange (from npp_nl.m)
dTf = (1/uf)*(ff*P*n - omg*Tf + 0.5*omg*(Tl + Te));
dTl = (1/uc)*((1-ff)*P*n + omg*Tf - 0.5*(2*M+omg)*Tl + 0.5*(2*M-omg)*Te);

dx = [dn; dc; dTf; dTl];
end
