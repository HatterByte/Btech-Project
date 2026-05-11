function dz = run_pid_only_nonlinear(t, z, p)
%RUN_PID_ONLY_NONLINEAR  Nonlinear PWR plant with PID control, no disturbance compensation.
%
%   Used to show the BASELINE performance (disturbance present, no IEID/ADRC).
%
%   States z (7x1):
%     z(1:4)  - [n; c; Tf; Tl] (plant core states)
%     z(5)    - rho_rod (rod reactivity)
%     z(6)    - x_int  (PID integrator)
%     z(7)    - x_d    (PID derivative filter)

% Unpack states
x_core  = z(1:4);
rho_rod = z(5);
x_int   = z(6);
x_d     = z(7);

n = x_core(1);

% Reference and error
if isfield(p, 'r_dev_fcn'); r_dev = p.r_dev_fcn(t); else; r_dev = 0.1 * (t >= 20); end
if isfield(p, 'n0'); n0 = p.n0; else; n0 = 0.5; end
y = n - n0;
e = r_dev - y;

% PID (no disturbance compensation — u_f is applied directly)
deriv_term = p.Kd * p.N * (e - x_d);
u_f = p.Kp * e + p.Ki * x_int + deriv_term;

% Disturbance acts on plant (but is NOT estimated or compensated)
if isfield(p, 'd_fcn'); d = p.d_fcn(t); else; d = disturbance_fcn(t); end

% Total reactivity rate input
u_total = u_f + d;

% Plant dynamics (nonlinear)
if isfield(p, 'Te_fcn'); Te = p.Te_fcn(t); else; Te = 290; end
if isfield(p, 'd_rho_fcn'); d_rho = p.d_rho_fcn(t); else; d_rho = 0; end
rho = total_reactivity(n, x_core(3), x_core(4), Te, rho_rod, d_rho);
dx_core  = pwr_nonlinear_dynamics(x_core, rho, Te);
drho_rod = p.Gr * u_total;

% PID states
dx_int = e;
dx_d   = -p.N * x_d + p.N * e;

dz = [dx_core; drho_rod; dx_int; dx_d];
end
