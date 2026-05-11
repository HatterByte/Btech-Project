function dz = run_adrc_closed_loop(t, z, p)
%RUN_ADRC_CLOSED_LOOP ODE wrapper for the ADRC controlled nonlinear PWR plant.
%
%   States z (8x1):
%     z(1:4) - [n; c; Tf; Tl] (Plant states)
%     z(5)   - rho_rod (Rod reactivity)
%     z(6:8) - [z1; z2; z3] (ESO states)

% Unpack plant states
x_plant = z(1:4);
rho_rod = z(5);
n = x_plant(1);

% Unpack ESO states
z_eso = z(6:8);

% Reference signal (deviation form)
if isfield(p, 'r_dev_fcn'); r_dev = p.r_dev_fcn(t); else; r_dev = 0.1 * (t >= 20); end

% Measured output (neutron density deviation)
if isfield(p, 'n0'); n0 = p.n0; else; n0 = 0.5; end
y = n - n0;

% External disturbance d(t) (reused from shared/disturbance_fcn.m)
if isfield(p, 'd_fcn'); d = p.d_fcn(t); else; d = disturbance_fcn(t); end

% Paper-matched ADRC control law.
u = adrc_control_law(r_dev, z_eso, p);

% Total Control Input to Plant = ADRC output + disturbance d
% Since d enters through the same channel as u (reactivity rate)
u_total = u + d;

% Plant Dynamics
if isfield(p, 'Te_fcn'); Te = p.Te_fcn(t); else; Te = 290; end
if isfield(p, 'd_rho_fcn'); d_rho = p.d_rho_fcn(t); else; d_rho = 0; end
rho = total_reactivity(n, x_plant(3), x_plant(4), Te, rho_rod, d_rho);
dx_plant = pwr_nonlinear_dynamics(x_plant, rho, Te);
drho_rod = p.Gr * u_total;

% ESO dynamics in the paper's reduced-order coordinates.
dz_eso = eso_dynamics(z_eso, y, u, p);

% Assemble derivative vector
dz = [dx_plant; drho_rod; dz_eso];

end
