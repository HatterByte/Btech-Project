function J = ieid_objective(gains)
%IEID_OBJECTIVE  Composite objective for IEID parameter tuning on linear PWR.
%
%   Designed specifically for the PWR context:
%     1. ITAE on disturbance window   — primary goal: reject d(t) after t=20s
%     2. Overshoot penalty            — safety: never exceed target power level
%     3. Steady-state error penalty   — converge accurately
%     4. Control effort penalty       — penalise saturation of rod velocity
%
%   INPUTS:
%     gains = [Kp, Ki, Kd, tau, log10(Q1)]   (5D)
%     Kd range [0, 0.2] — allows 0 if optimizer chooses it
%
%   OUTPUT:
%     J : scalar cost (lower = better)

Kp  = gains(1);
Ki  = gains(2);
Kd  = gains(3);
tau = gains(4);
Q1  = 10^gains(5);

%% Build param struct
p = ieid_opt_params(Kp, Ki, Kd, tau, Q1);
if isempty(p); J = 1e6; return; end

%% Simulate
tspan = [0 60];
z0    = zeros(14, 1);
opts  = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'MaxStep', 0.1);

try
    [t, z] = ode15s(@(t,z) ieid_opt_ode(t, z, p), tspan, z0, opts);
catch
    J = 1e6; return;
end

%% Signals
y = (p.C * z(:, 1:5)')';   % neutron density deviation
r = 0.1 * (t >= 20);       % reference
e = r - y;                  % tracking error

% Control effort reconstruction (for saturation penalty)
u_sig = zeros(size(t));
for k = 1:numel(t)
    ek    = e(k);
    xi    = z(k, 11);
    xd    = z(k, 12);
    xf1   = z(k, 13);
    u_f   = Kp*ek + Ki*xi + Kd*100*(ek - xd);   % N=100 derivative filter
    u_sig(k) = u_f - xf1;
end

% Divergence guard
if any(~isfinite(y)) || max(abs(y)) > 10
    J = 1e6; return;
end

%% --- Component 1: ITAE on DISTURBANCE REJECTION window (t >= 20s) ---
% This is the heart of IEID's purpose — reject d(t) while tracking r.
% Weighted by time to penalise slow recovery more than fast initial error.
mask_dist = t >= 20;
t_dist = t(mask_dist) - 20;   % reset time origin to step onset
e_dist = e(mask_dist);
J_dist = trapz(t_dist, t_dist .* abs(e_dist));

%% --- Component 2: Overshoot penalty ---
% Nuclear safety: output must NEVER exceed reference by more than 0.5%.
target   = 0.1;
overshoot = max(0, max(y) - target);   % positive when output exceeds target
J_over   = 500 * overshoot^2;          % quadratic, high weight

%% --- Component 3: Steady-state error (last 5 seconds) ---
ss_mask = t >= 55;
J_ss    = 200 * mean(abs(e(ss_mask)));

%% --- Component 4: Control effort saturation penalty ---
% Rod velocity limit: ±0.01 δk/k/s. Penalise time spent at saturation.
u_lim = 0.01;
J_sat = 5 * trapz(t, max(0, abs(u_sig) - u_lim).^2);

%% --- Composite ---
% Weights reflect priority: disturbance rejection > SS accuracy > safety > sat
J = J_dist ...
  + J_over ...
  + J_ss   ...
  + J_sat;

end
