function dz = ieid_opt_ode(t, z, p)
%IEID_OPT_ODE  Compact IEID closed-loop ODE for use inside optimizer.
%   Identical physics to run_closed_loop / run_ieid_closed_loop but
%   operates on the LINEAR plant (no disturbance in nominal case,
%   disturbance_fcn active in disturbance case).
%
%   States z (14x1):
%     z(1:5)   plant deviation states [dn; dc; dTf; dTl; drho_r]
%     z(6:10)  observer state x_hat
%     z(11)    PID integrator x_int
%     z(12)    PID derivative filter x_d
%     z(13:14) IEID LPF states [xf1; xf2]

x     = z(1:5);
x_hat = z(6:10);
x_int = z(11);
x_d   = z(12);
xf1   = z(13);
xf2   = z(14);

% Output deviation
y = p.C * x;

% Reference
r = 0.1 * (t >= 20);
e = r - y;

% PID
N   = p.N;
u_f = p.Kp*e + p.Ki*x_int + p.Kd*N*(e - x_d);

% IEID compensation
d_tilde = xf1;
u = u_f - d_tilde;

% Disturbance (active during optimization for disturbance rejection tuning)
d = disturbance_fcn(t);

% Plant
dx = p.A_real * x + p.B_real * u + p.B_d * d;

% Observer
dx_hat = p.A_obs * x_hat + p.B_obs * [u_f; y];

% IEID raw estimate: standard form d_hat = Bnp * L * (y - y_hat) + (u_f - u)
y_hat = p.C * x_hat;
d_hat = p.Bnp * (p.L_obs * (y - y_hat)) + (u_f - u);

% LPF
tau   = p.tau;
dxf1  = xf2;
dxf2  = (d_hat - xf1 - 1.41*tau*xf2) / (tau^2);

% PID filter states
dx_int = e;
dx_d   = -N*x_d + N*e;

dz = [dx; dx_hat; dx_int; dx_d; dxf1; dxf2];
end
