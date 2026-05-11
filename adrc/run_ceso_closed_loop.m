function dz = run_ceso_closed_loop(t, z, p)
%RUN_CESO_CLOSED_LOOP ODE wrapper for CESO-ADRC controlled nonlinear PWR.
%
%   States z (11x1):
%     z(1:4)  - [n; c; Tf; Tl] plant states
%     z(5)    - rho_rod
%     z(6:11) - two 3-state CESO levels

x_plant = z(1:4);
rho_rod = z(5);
n = x_plant(1);
z_ceso = z(6:11);

if isfield(p, 'r_dev_fcn'); r_dev = p.r_dev_fcn(t); else; r_dev = 0.1 * (t >= 20); end
if isfield(p, 'n0'); n0 = p.n0; else; n0 = 0.5; end
y = n - n0;
if isfield(p, 'd_fcn'); d = p.d_fcn(t); else; d = disturbance_fcn(t); end

z1 = z_ceso(1:3);
z2 = z_ceso(4:6);
x_hat = [z2(1); z2(2); z2(3) + z1(3)];
u = adrc_control_law(r_dev, x_hat, p);
u_total = u + d;

if isfield(p, 'Te_fcn'); Te = p.Te_fcn(t); else; Te = 290; end
if isfield(p, 'd_rho_fcn'); d_rho = p.d_rho_fcn(t); else; d_rho = 0; end
rho = total_reactivity(n, x_plant(3), x_plant(4), Te, rho_rod, d_rho);
dx_plant = pwr_nonlinear_dynamics(x_plant, rho, Te);
drho_rod = p.Gr * u_total;

dz_ceso = ceso_dynamics(z_ceso, y, u, p);

dz = [dx_plant; drho_rod; dz_ceso];
end
