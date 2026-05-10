%% ================================================================
%%  MAIN_IEID_PWR  —  IEID-Based Disturbance Rejection for PWR
%% ================================================================
%
%  This script replicates the Simulink closed-loop simulation using
%  pure MATLAB (no Simulink required).
%
%  System:
%    - PWR 5-state linearised plant (A_n, B_n, C)
%    - PID outer-loop controller (Kp, Ki, Kd, N=100)
%    - Full-state Luenberger observer
%    - IEID estimator + 2nd-order low-pass filter
%
%  Augmented ODE state  z (14×1):
%    z(1:5)   plant states  x
%    z(6:10)  observer states  x_hat
%    z(11)    PID integral state
%    z(12)    PID derivative filter state  (N = 100)
%    z(13)    IEID LPF state 1  (= d_tilde, filter output)
%    z(14)    IEID LPF state 2  (filter rate)
%
%  Reference: step from 0 to 0.1 (deviation) at t = 20 s.
%  Disturbance: d(t) = 0.002*(sin(5t) + 0.5*tanh(4(t-20)) + cos(3t))
%
%  Run:  >> main_ieid_pwr
%
%% ================================================================

clear; clc; close all;

%% ---- Load system parameters (params_pwr.m) ----------------------
params_pwr;     % defines A_n, B_n, C, D, A_real, B_real, B_d,
                %         L, Bnp, tau, numF, denF, Kp, Ki, Kd,
                %         A_obs, B_obs, C_obs, D_obs

%% ---- Verify observer stability (sanity check) -------------------
obs_poles = eig(A_n - L * C);
fprintf('\n--- Observer poles (should all be stable) ---\n');
for k = 1:length(obs_poles)
    fprintf('  pole %d : %.4f + %.4fi\n', k, real(obs_poles(k)), imag(obs_poles(k)));
end
if any(real(obs_poles) > 0)
    warning('Observer is UNSTABLE — check L gain!');
else
    fprintf('Observer is stable. ✓\n\n');
end

%% ---- PID derivative filter coefficient --------------------------
N_pid = 100;    % Simulink default derivative filter coefficient

%% ---- Pack parameters into struct --------------------------------
p.A_real = A_real;
p.B_real = B_real;
p.B_d    = B_d;
p.A_obs  = A_obs;          % = A_n - L*C
p.B_obs  = B_obs;          % = [B_n, L]
p.C      = C;
p.L      = L;
p.Bnp    = Bnp;
p.Kp     = Kp;
p.Ki     = Ki;
p.Kd     = Kd;
p.N      = N_pid;
p.tau    = tau;

%% ---- Initial conditions -----------------------------------------
%  All states start at zero (nominal operating point, deviation form)
z0 = zeros(14, 1);

%% ---- Simulation time span ---------------------------------------
t_span = [0, 60];          % 60 seconds; step occurs at t = 20 s

%% ---- ODE solver options -----------------------------------------
%  ode15s is preferred: the PWR system + observer can be stiff
opts = odeset('RelTol', 1e-7, ...
              'AbsTol', 1e-9 * ones(14,1), ...
              'MaxStep', 0.05);

%% ---- Solve closed-loop ODE -------------------------------------
fprintf('Running ode15s simulation ...\n');
tic;
[t_sol, Z_sol] = ode15s(@(t,z) run_closed_loop(t, z, p), t_span, z0, opts);
elapsed = toc;
fprintf('Simulation complete in %.2f s  (%d time points).\n\n', elapsed, length(t_sol));

%% ---- Post-process: recompute all signals at solution times ------
N_pts = length(t_sol);

sig.r       = zeros(N_pts, 1);
sig.r_dev   = zeros(N_pts, 1);
sig.y       = zeros(N_pts, 1);
sig.e       = zeros(N_pts, 1);
sig.u_f     = zeros(N_pts, 1);
sig.u       = zeros(N_pts, 1);
sig.d       = zeros(N_pts, 1);
sig.d_hat   = zeros(N_pts, 1);
sig.d_tilde = zeros(N_pts, 1);

for k = 1:N_pts
    t_k = t_sol(k);
    z_k = Z_sol(k, :)';

    % Unpack state
    x_plant = z_k(1:5);
    x_hat   = z_k(6:10);
    x_int   = z_k(11);   %#ok<NASGU>
    x_d     = z_k(12);
    xf1     = z_k(13);   % = d_tilde

    % Reference
    r_k     = 0.5 + 0.1 * (t_k >= 20);
    r_dev_k = r_k - 0.5;

    % Plant output
    y_k = C * x_plant;

    % Error
    e_k = r_dev_k - y_k;

    % PID
    deriv_k = Kd * N_pid * (e_k - x_d);
    u_f_k   = Kp * e_k + Ki * x_int + deriv_k;

    % Compensation
    d_tilde_k = xf1;
    u_k       = u_f_k - d_tilde_k;

    % Disturbance
    d_k = disturbance_fcn(t_k);

    % IEID raw estimate
    d_hat_k = ieid_estimator_fcn(u_f_k, u_k, x_hat, y_k, L, C, Bnp);

    % Store
    sig.r(k)       = r_k;
    sig.r_dev(k)   = r_dev_k;
    sig.y(k)       = y_k;
    sig.e(k)       = e_k;
    sig.u_f(k)     = u_f_k;
    sig.u(k)       = u_k;
    sig.d(k)       = d_k;
    sig.d_hat(k)   = d_hat_k;
    sig.d_tilde(k) = d_tilde_k;
end

%% ---- Print steady-state summary ---------------------------------
fprintf('--- Steady-state summary (t = 55–60 s) ---\n');
ss_idx = t_sol >= 55;
fprintf('  y     mean = %.6f  (target r_dev = 0.1)\n', mean(sig.y(ss_idx)));
fprintf('  e     mean = %.2e\n', mean(sig.e(ss_idx)));
fprintf('  d     mean = %.6f\n', mean(sig.d(ss_idx)));
fprintf('  d̃     mean = %.6f\n', mean(sig.d_tilde(ss_idx)));
fprintf('  |d - d̃| mean = %.2e\n', mean(abs(sig.d(ss_idx) - sig.d_tilde(ss_idx))));

%% ---- Plot results -----------------------------------------------
plot_ieid_results(t_sol, sig);

fprintf('\nDone. Four figure windows generated.\n');
