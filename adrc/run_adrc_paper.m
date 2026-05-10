function dz = run_adrc_paper(t, z, p)
%RUN_ADRC_PAPER  Standard 3rd-order ESO ADRC on the paper's nonlinear PWR plant.
%
%   Uses corrected plant parameters confirmed from paper's Simulink source:
%     - Tf0=650, Tl0=314, Te0=290  (from tot_react.m)
%     - af=(n-4.24)*1e-5           (from tot_react.m)
%     - B1=0.8*B, lmb1=0.8*lmb    (from npp_nl.m, uncertain model)
%
%   Architecture: standard ADRC with 3rd-order LESO.
%   The ESO observes n_r as output and estimates [n_r, dn_r/dt, total_disturbance].
%
%   States z (8x1):
%     z(1:4) = [n; c; Tf; Tl]  — plant core states
%     z(5)   = rho_rod          — control rod reactivity
%     z(6:8) = [z1; z2; z3]    — ESO states (3rd-order)

%% Unpack states
n       = z(1);
rho_rod = z(5);
z_eso   = z(6:8);
z1_hat  = z_eso(1);
z2_hat  = z_eso(2);
z3_hat  = z_eso(3);

%% Reference signal — Scenario 1: power tracking
if     t < 50;   n_ref = 1.0;
elseif t < 100;  n_ref = 0.9;
elseif t < 150;  n_ref = 0.7;
elseif t < 156;  n_ref = 0.7 + 0.05*(t - 150);
else;            n_ref = 1.0;
end

%% Coolant inlet temperature — Scenario 2
Te = p.Te0;
if t >= 200; Te = Te - 2.0; end
if t >= 240; Te = Te + min(2.0*(t - 240), 7.0); end

%% External reactivity disturbance — Scenario 3
d_rho = 0;
if t >= 290; d_rho = -0.0002; end
if t >= 320; d_rho = d_rho + min(0.0002*(t - 320), 0.002); end

%% Measured output (n_r is the plant output)
y = n;

%% ADRC control law (standard, in original output space)
u_bar = p.k1*(n_ref - z1_hat) - p.k2*z2_hat;
u_raw = (u_bar - z3_hat) / p.b0;
u     = max(p.u_min, min(p.u_max, u_raw));

%% Plant dynamics
dx_core  = pwr_plant_paper(z(1:4), rho_rod, Te, d_rho, p.Tf_ref, p.Tl_ref, p.Te_ref);
drho_rod = p.Gr * u;

%% ESO dynamics (3rd-order LESO)
dz_eso = eso_dynamics(z_eso, y, u, p.wo, p.b0);

%% Assemble
dz = [dx_core; drho_rod; dz_eso];
end
