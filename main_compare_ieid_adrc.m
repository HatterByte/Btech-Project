%% MAIN_COMPARE_IEID_ADRC
% Comparative study of IEID and ADRC on a nonlinear PWR plant.
% Includes PID-only baseline (with disturbance, no compensation) for
% BOTH the old linearized plant and the new nonlinear plant.

clear; clc; close all;

% Add folders to path
addpath('plant');
addpath('ieid');
addpath('adrc');
addpath('shared');

% ODE options (same for all simulations)
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
tspan = [0 60];

%% =========================================================================
%% SECTION 1 — NEW NONLINEAR PLANT
%% =========================================================================
p_nl = get_nonlinear_params();

% Common initial conditions for nonlinear plant
% [n; c; Tf; Tl; rho_rod] at equilibrium
x0_nl = [0.5; 0.5; 650; 314; 0];

%% --- 1a. PID-Only Baseline (nonlinear plant, no compensation) ---
fprintf('Running PID-Only baseline (nonlinear plant)...\n');
z0_base_nl = [x0_nl; 0; 0];   % [plant(5); x_int; x_d]
[t_base_nl, z_base_nl] = ode15s(@(t,z) run_pid_only_nonlinear(t, z, p_nl), ...
                                  tspan, z0_base_nl, opts);

res_base_nl.t = t_base_nl;
res_base_nl.y = z_base_nl(:,1) - 0.5;
res_base_nl.r = 0.1 * (t_base_nl >= 20);
res_base_nl.u = zeros(length(t_base_nl), 1);
for k = 1:length(t_base_nl)
    y   = z_base_nl(k,1) - 0.5;
    e   = res_base_nl.r(k) - y;
    xi  = z_base_nl(k,6);
    xd  = z_base_nl(k,7);
    res_base_nl.u(k) = p_nl.Kp*e + p_nl.Ki*xi + p_nl.Kd*p_nl.N*(e - xd);
end

%% --- 1b. IEID on nonlinear plant ---
fprintf('Running IEID simulation (nonlinear plant)...\n');
z0_ieid = [x0_nl; zeros(5,1); 0; 0; 0; 0];   % [plant(5); obs(5); x_int; x_d; lpf(2)]
[t_ieid, z_ieid] = ode15s(@(t,z) run_ieid_closed_loop(t, z, p_nl), ...
                            tspan, z0_ieid, opts);

res_ieid.t     = t_ieid;
res_ieid.y     = z_ieid(:,1) - 0.5;
res_ieid.r     = 0.1 * (t_ieid >= 20);
res_ieid.d     = arrayfun(@disturbance_fcn, t_ieid);
res_ieid.d_est = z_ieid(:,13);   % d_tilde (LPF output)
res_ieid.u     = zeros(length(t_ieid), 1);
for k = 1:length(t_ieid)
    y       = z_ieid(k,1) - 0.5;
    e       = res_ieid.r(k) - y;
    xi      = z_ieid(k,11);
    xd      = z_ieid(k,12);
    dtilde  = z_ieid(k,13);
    u_f     = p_nl.Kp*e + p_nl.Ki*xi + p_nl.Kd*p_nl.N*(e - xd);
    res_ieid.u(k) = u_f - dtilde;
end

%% --- 1c. ADRC on nonlinear plant ---
fprintf('Running ADRC simulation (nonlinear plant)...\n');
z0_adrc = [x0_nl; zeros(3,1)];   % [plant(5); eso(3)]
[t_adrc, z_adrc] = ode15s(@(t,z) run_adrc_closed_loop(t, z, p_nl), ...
                            tspan, z0_adrc, opts);

res_adrc.t     = t_adrc;
res_adrc.y     = z_adrc(:,1) - 0.5;
res_adrc.r     = 0.1 * (t_adrc >= 20);
res_adrc.d     = arrayfun(@disturbance_fcn, t_adrc);
res_adrc.d_est = z_adrc(:,8) / p_nl.b0;   % z3 normalised to input space
res_adrc.u     = zeros(length(t_adrc), 1);
for k = 1:length(t_adrc)
    z_eso = z_adrc(k,6:8)';
    r_dev = 0.1 * (t_adrc(k) >= 20);
    u_bar = p_nl.k1*(r_dev - z_eso(1)) - p_nl.k2*z_eso(2);
    res_adrc.u(k) = (u_bar - z_eso(3)) / p_nl.b0;
end

%% =========================================================================
%% SECTION 2 — OLD LINEARIZED PLANT
%% =========================================================================
p_lin = get_linear_params();

% Common initial conditions: deviation form, starts at equilibrium (all zeros)
x0_lin = zeros(5, 1);

%% --- 2a. PID-Only Baseline (old linear plant, no compensation) ---
fprintf('Running PID-Only baseline (old linear plant)...\n');
z0_base_lin = [x0_lin; 0; 0];   % [plant(5); x_int; x_d]
[t_base_lin, z_base_lin] = ode15s(@(t,z) run_pid_only_linear(t, z, p_lin), ...
                                    tspan, z0_base_lin, opts);

res_base_lin.t = t_base_lin;
res_base_lin.y = (p_lin.C * z_base_lin(:, 1:5)')';   % neutron density deviation
res_base_lin.r = 0.1 * (t_base_lin >= 20);
res_base_lin.u = zeros(length(t_base_lin), 1);
for k = 1:length(t_base_lin)
    y  = res_base_lin.y(k);
    e  = res_base_lin.r(k) - y;
    xi = z_base_lin(k,6);
    xd = z_base_lin(k,7);
    res_base_lin.u(k) = p_lin.Kp*e + p_lin.Ki*xi + p_lin.Kd*p_lin.N*(e - xd);
end

%% --- 2b. IEID on old linear plant ---
fprintf('Running IEID simulation (old linear plant)...\n');
z0_ieid_lin = zeros(14, 1);   % [plant(5); obs(5); x_int; x_d; lpf(2)]

% Build old-plant param struct for run_closed_loop (old ODE function)
p_lin_ieid        = p_lin;
p_lin_ieid.A_real = p_lin.A_real;
p_lin_ieid.B_real = p_lin.B_real;
p_lin_ieid.B_d    = p_lin.B_d;
p_lin_ieid.A_obs  = p_lin.A_obs;
p_lin_ieid.B_obs  = p_lin.B_obs;
p_lin_ieid.C      = p_lin.C;
p_lin_ieid.L      = p_lin.L_obs;
p_lin_ieid.Bnp    = p_lin.Bnp;
p_lin_ieid.tau    = p_lin.tau;

[t_ieid_lin, z_ieid_lin] = ode15s(@(t,z) run_closed_loop(t, z, p_lin_ieid), ...
                                    tspan, z0_ieid_lin, opts);

res_ieid_lin.t     = t_ieid_lin;
res_ieid_lin.y     = (p_lin.C * z_ieid_lin(:, 1:5)')';
res_ieid_lin.r     = 0.1 * (t_ieid_lin >= 20);
res_ieid_lin.d     = arrayfun(@disturbance_fcn, t_ieid_lin);
res_ieid_lin.d_est = z_ieid_lin(:,13);
res_ieid_lin.u     = zeros(length(t_ieid_lin), 1);
for k = 1:length(t_ieid_lin)
    y      = res_ieid_lin.y(k);
    e      = res_ieid_lin.r(k) - y;
    xi     = z_ieid_lin(k,11);
    xd     = z_ieid_lin(k,12);
    dtilde = z_ieid_lin(k,13);
    u_f    = p_lin.Kp*e + p_lin.Ki*xi + p_lin.Kd*p_lin.N*(e - xd);
    res_ieid_lin.u(k) = u_f - dtilde;
end

%% =========================================================================
%% SECTION 3 — METRICS
%% =========================================================================
m_ieid     = metrics(t_ieid,     res_ieid.y,     res_ieid.r,     res_ieid.u);
m_adrc     = metrics(t_adrc,     res_adrc.y,     res_adrc.r,     res_adrc.u);
m_base_nl  = metrics(t_base_nl,  res_base_nl.y,  res_base_nl.r,  res_base_nl.u);
m_base_lin = metrics(t_base_lin, res_base_lin.y, res_base_lin.r, res_base_lin.u);
m_ieid_lin = metrics(t_ieid_lin, res_ieid_lin.y, res_ieid_lin.r, res_ieid_lin.u);

fprintf('\n=== Performance Comparison (Nonlinear Plant) ===\n');
fprintf('%-22s | %-12s | %-12s | %-12s\n', 'Metric', 'PID-Only', 'IEID', 'ADRC');
fprintf('%s\n', repmat('-',1,65));
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'RMS Error',       m_base_nl.rms_error,   m_ieid.rms_error,   m_adrc.rms_error);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Rise Time (s)',    m_base_nl.rise_time,   m_ieid.rise_time,   m_adrc.rise_time);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Settling Time (s)',m_base_nl.settling_time,m_ieid.settling_time,m_adrc.settling_time);
fprintf('%-22s | %-11.2f%% | %-11.2f%% | %-11.2f%%\n', 'Overshoot',  m_base_nl.overshoot,   m_ieid.overshoot,   m_adrc.overshoot);
fprintf('%-22s | %-12.4e | %-12.4e | %-12.4e\n', 'SS Error',         m_base_nl.sse,         m_ieid.sse,         m_adrc.sse);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Control Energy',   m_base_nl.control_energy,m_ieid.control_energy,m_adrc.control_energy);

fprintf('\n=== Performance Comparison (Old Linear Plant) ===\n');
fprintf('%-22s | %-12s | %-12s\n', 'Metric', 'PID-Only', 'IEID');
fprintf('%s\n', repmat('-',1,50));
fprintf('%-22s | %-12.4f | %-12.4f\n', 'RMS Error',        m_base_lin.rms_error,    m_ieid_lin.rms_error);
fprintf('%-22s | %-12.4f | %-12.4f\n', 'Rise Time (s)',     m_base_lin.rise_time,    m_ieid_lin.rise_time);
fprintf('%-22s | %-12.4f | %-12.4f\n', 'Settling Time (s)', m_base_lin.settling_time,m_ieid_lin.settling_time);
fprintf('%-22s | %-11.2f%% | %-11.2f%%\n', 'Overshoot',     m_base_lin.overshoot,   m_ieid_lin.overshoot);
fprintf('%-22s | %-12.4e | %-12.4e\n', 'SS Error',          m_base_lin.sse,          m_ieid_lin.sse);
fprintf('%-22s | %-12.4f | %-12.4f\n', 'Control Energy',    m_base_lin.control_energy,m_ieid_lin.control_energy);

%% =========================================================================
%% SECTION 4 — PLOTS
%% =========================================================================
plotting(res_ieid, res_adrc, res_base_nl, res_ieid_lin, res_base_lin);
