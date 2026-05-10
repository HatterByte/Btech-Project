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
z1_hat = z_eso(1);
z2_hat = z_eso(2);
z3_hat = z_eso(3);

% Reference signal (deviation form)
r_step = 0.5 + 0.1 * (t >= 20);
r_dev  = r_step - 0.5;

% Measured output (neutron density deviation)
y = n - 0.5;

% External disturbance d(t) (reused from shared/disturbance_fcn.m)
d = disturbance_fcn(t);

% ADRC Control Law
% u_bar = k1*(r_dev - z1_hat) - k2*z2_hat
u_bar = p.k1 * (r_dev - z1_hat) - p.k2 * z2_hat;
% u = (u_bar - z3_hat) / b0
u = (u_bar - z3_hat) / p.b0;

% Total Control Input to Plant = ADRC output + disturbance d
% Since d enters through the same channel as u (reactivity rate)
u_total = u + d;

% Plant Dynamics
Te = 290; % Constant inlet temperature
rho = total_reactivity(n, x_plant(3), x_plant(4), Te, rho_rod);
dx_plant = pwr_nonlinear_dynamics(x_plant, rho, Te);
drho_rod = p.Gr * u_total;

% ESO Dynamics
dz_eso = eso_dynamics(z_eso, y, u, p.wo, p.b0);

% Assemble derivative vector
dz = [dx_plant; drho_rod; dz_eso];

end
