function dz = run_ieid_closed_loop(t, z, p)
%RUN_IEID_CLOSED_LOOP ODE wrapper for IEID controlled nonlinear PWR plant.
%
%   States z (14x1):
%     z(1:4)   - [n; c; Tf; Tl] (Plant states)
%     z(5)     - rho_rod (Rod reactivity)
%     z(6:10)  - x_hat (Linear observer states)
%     z(11)    - x_int (PID integrator)
%     z(12)    - x_d (PID derivative filter)
%     z(13:14) - [xf1; xf2] (IEID LPF states)

% Unpack plant states
x_plant_core = z(1:4);
rho_rod = z(5);
n = x_plant_core(1);

% Unpack observer/controller states
x_hat = z(6:10);
x_int = z(11);
x_d   = z(12);
xf1   = z(13); % d_tilde
xf2   = z(14);

% Reference
r_step = 0.5 + 0.1 * (t >= 20);
r_dev  = r_step - 0.5;

% Output deviation
y = n - 0.5;

% Tracking error
e = r_dev - y;

% PID Controller (derivative filter coefficient N from params)
N = p.N;
deriv_term = p.Kd * N * (e - x_d);
u_f = p.Kp * e + p.Ki * x_int + deriv_term;

% IEID Compensation
d_tilde = xf1;
u = u_f - d_tilde;

% External disturbance
d = disturbance_fcn(t);

% Total reactivity rate input to plant
u_total = u + d;

% --- Plant Dynamics (Nonlinear) ---
Te = 290;
rho = total_reactivity(n, x_plant_core(3), x_plant_core(4), Te, rho_rod);
dx_plant_core = pwr_nonlinear_dynamics(x_plant_core, rho, Te);
drho_rod = p.Gr * u_total;

% --- Controller/Observer Dynamics ---

% 1. IEID Estimator (Raw estimate d_hat)
d_hat = ieid_estimator_fcn(u_f, u, x_hat, y, p.L_ieid, p.C, pinv(p.B));


% 2. Observer: dx_hat = (A - L*C)*x_hat + B*u_f + L*y
A_obs = p.A - p.L_ieid * p.C;
dx_hat = A_obs * x_hat + p.B * u_f + p.L_ieid * y;

% 3. PID integrator and derivative filter
dx_int = e;
dx_d   = -p.N * x_d + p.N * e;

% 4. IEID 2nd-order LPF
tau = 0.02;
dxf1 = xf2;
dxf2 = (d_hat - xf1 - 1.41 * tau * xf2) / (tau^2);

% Assemble derivative vector
dz = [dx_plant_core; drho_rod; dx_hat; dx_int; dx_d; dxf1; dxf2];

end
