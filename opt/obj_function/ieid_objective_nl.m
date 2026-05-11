function J = ieid_objective_nl(gains)
%IEID_OBJECTIVE_NL Composite objective for IEID parameter tuning on NONLINEAR PWR.
%
%   INPUTS:
%     gains = [Kp, Ki, Kd, tau, log10(Q1)]   (5D)
%
%   OUTPUT:
%     J : scalar cost (lower = better)

Kp  = gains(1);
Ki  = gains(2);
Kd  = gains(3);
tau = gains(4);
Q1  = 10^gains(5);

%% Build param struct from nonlinear defaults
p = get_nonlinear_params();
p.Kp = Kp;
p.Ki = Ki;
p.Kd = Kd;
p.tau = tau;

% Recompute observer gain with new Q
Q = diag([Q1, 1e2, 10, 10, 1]);
try
    p.L_ieid = lqe(p.A, eye(5), p.C, Q, 1);
catch
    J = 1e6; return;
end

% Set up disturbance and reference functions for the simulation
p.r_dev_fcn = @(t) 0.1 * (t >= 20);
p.d_fcn     = @disturbance_fcn;

%% Simulate
tspan = [0 60];
x0_nl = [0.5; 0.5; 650; 314; 0];
z0    = [x0_nl; zeros(5,1); 0; 0; 0; 0];
% Looser tolerances than the linear case to speed up GCRA over nonlinear ODE
opts  = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'MaxStep', 0.2);

try
    [t, z] = ode15s(@(t,z) run_ieid_closed_loop(t, z, p), tspan, z0, opts);
catch
    J = 1e6; return;
end

%% Signals
y = z(:, 1) - 0.5;         % neutron density deviation
r = arrayfun(p.r_dev_fcn, t);
e = r - y;

% Divergence guard
if any(~isfinite(y)) || max(abs(y)) > 2
    J = 1e6; return;
end

% Control effort reconstruction (for saturation penalty)
u_sig = zeros(size(t));
for k = 1:numel(t)
    ek    = e(k);
    xi    = z(k, 11);
    xd    = z(k, 12);
    xf1   = z(k, 13);
    u_f   = p.Kp*ek + p.Ki*xi + p.Kd*p.N*(ek - xd);
    u_sig(k) = u_f - xf1;
end

%% --- Component 1: ITAE on DISTURBANCE REJECTION window (t >= 20s) ---
mask_dist = t >= 20;
t_dist = t(mask_dist) - 20;   
e_dist = e(mask_dist);
J_dist = trapz(t_dist, t_dist .* abs(e_dist));

%% --- Component 2: Overshoot penalty ---
% Positive when output exceeds 0.1
overshoot = max(0, max(y) - 0.1);   
J_over   = 500 * overshoot^2;       

%% --- Component 3: Steady-state error (last 5 seconds) ---
ss_mask = t >= 55;
J_ss    = 200 * mean(abs(e(ss_mask)));

%% --- Component 4: Control effort saturation penalty ---
u_lim = p.u_max; % 0.01
J_sat = 5 * trapz(t, max(0, abs(u_sig) - u_lim).^2);

%% --- Composite ---
J = J_dist + J_over + J_ss + J_sat;

end
