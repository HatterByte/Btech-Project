function dz = run_pid_only_linear(t, z, p)
%RUN_PID_ONLY_LINEAR  Old linearized PWR plant with PID control, no disturbance compensation.
%
%   Used to show the BASELINE performance (disturbance present, no IEID).
%
%   States z (7x1):
%     z(1:5)  - plant deviation states [dn; dc; dTf; dTl; drho_r]
%     z(6)    - x_int  (PID integrator)
%     z(7)    - x_d    (PID derivative filter)

% Unpack
x_plant = z(1:5);
x_int   = z(6);
x_d     = z(7);

% Output (neutron density deviation)
y = p.C * x_plant;

% Reference deviation and error
r_dev = 0.1 * (t >= 20);
e     = r_dev - y;

% PID (no disturbance compensation)
deriv_term = p.Kd * p.N * (e - x_d);
u_f = p.Kp * e + p.Ki * x_int + deriv_term;

% Disturbance (NOT compensated)
d = disturbance_fcn(t);

% Plant: ẋ = A_real*x + B_real*u_f + B_d*d
dx_plant = p.A_real * x_plant + p.B_real * u_f + p.B_d * d;

% PID states
dx_int = e;
dx_d   = -p.N * x_d + p.N * e;

dz = [dx_plant; dx_int; dx_d];
end
