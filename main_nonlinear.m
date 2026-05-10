%% MAIN_NONLINEAR
% Simulation of the NEW NONLINEAR PWR plant.
% Compares:
%   1. PID-only (with disturbance, no compensation) — baseline
%   2. IEID (disturbance estimation + rejection)
%   3. ADRC (ESO-based active disturbance rejection)

clear; clc; close all;

addpath('plant');
addpath('ieid');
addpath('adrc');
addpath('shared');

%% Parameters
p = get_nonlinear_params();

tspan = [0 60];
opts  = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Common initial conditions: equilibrium point [n; c; Tf; Tl; rho_rod]
x0 = [0.5; 0.5; 650; 314; 0];

%% =========================================================================
%% 1. PID-Only Baseline (disturbance present, no compensation)
%% =========================================================================
fprintf('Running PID-Only baseline (nonlinear plant)...\n');
z0_base = [x0; 0; 0];     % [plant(5); x_int; x_d]

[t_base, z_base] = ode15s(@(t,z) run_pid_only_nonlinear(t, z, p), tspan, z0_base, opts);

res_base.t = t_base;
res_base.y = z_base(:, 1) - 0.5;    % neutron density deviation
res_base.r = 0.1 * (t_base >= 20);
res_base.u = zeros(length(t_base), 1);
for k = 1:length(t_base)
    y  = res_base.y(k);
    e  = res_base.r(k) - y;
    xi = z_base(k, 6);
    xd = z_base(k, 7);
    res_base.u(k) = p.Kp*e + p.Ki*xi + p.Kd*p.N*(e - xd);
end

%% =========================================================================
%% 2. IEID (disturbance estimation and rejection)
%% =========================================================================
fprintf('Running IEID (nonlinear plant)...\n');
z0_ieid = [x0; zeros(5,1); 0; 0; 0; 0];   % [plant(5); obs(5); x_int; x_d; lpf(2)]

[t_ieid, z_ieid] = ode15s(@(t,z) run_ieid_closed_loop(t, z, p), tspan, z0_ieid, opts);

res_ieid.t     = t_ieid;
res_ieid.y     = z_ieid(:, 1) - 0.5;
res_ieid.r     = 0.1 * (t_ieid >= 20);
res_ieid.d     = arrayfun(@disturbance_fcn, t_ieid);
res_ieid.d_est = z_ieid(:, 13);    % d_tilde (IEID LPF output)
res_ieid.u     = zeros(length(t_ieid), 1);
for k = 1:length(t_ieid)
    y      = res_ieid.y(k);
    e      = res_ieid.r(k) - y;
    xi     = z_ieid(k, 11);
    xd     = z_ieid(k, 12);
    dtilde = z_ieid(k, 13);
    u_f    = p.Kp*e + p.Ki*xi + p.Kd*p.N*(e - xd);
    res_ieid.u(k) = u_f - dtilde;
end

%% =========================================================================
%% 3. ADRC (ESO-based active disturbance rejection)
%% =========================================================================
fprintf('Running ADRC (nonlinear plant)...\n');
z0_adrc = [x0; zeros(3,1)];   % [plant(5); eso(3)]

[t_adrc, z_adrc] = ode15s(@(t,z) run_adrc_closed_loop(t, z, p), tspan, z0_adrc, opts);

res_adrc.t     = t_adrc;
res_adrc.y     = z_adrc(:, 1) - 0.5;
res_adrc.r     = 0.1 * (t_adrc >= 20);
res_adrc.d     = arrayfun(@disturbance_fcn, t_adrc);
res_adrc.d_est = z_adrc(:, 8) / p.b0;    % z3 normalised to input space
res_adrc.u     = zeros(length(t_adrc), 1);
for k = 1:length(t_adrc)
    z_eso = z_adrc(k, 6:8)';
    r_dev = 0.1 * (t_adrc(k) >= 20);
    u_bar = p.k1*(r_dev - z_eso(1)) - p.k2*z_eso(2);
    res_adrc.u(k) = (u_bar - z_eso(3)) / p.b0;
end

%% =========================================================================
%% 4. Metrics
%% =========================================================================
m_base = metrics(t_base, res_base.y, res_base.r, res_base.u);
m_ieid = metrics(t_ieid, res_ieid.y, res_ieid.r, res_ieid.u);
m_adrc = metrics(t_adrc, res_adrc.y, res_adrc.r, res_adrc.u);

fprintf('\n=== Performance Comparison — Nonlinear PWR Plant ===\n');
fprintf('%-22s | %-12s | %-12s | %-12s\n', 'Metric', 'PID-Only', 'IEID', 'ADRC');
fprintf('%s\n', repmat('-', 1, 65));
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'RMS Error',        m_base.rms_error,    m_ieid.rms_error,    m_adrc.rms_error);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Rise Time (s)',     m_base.rise_time,    m_ieid.rise_time,    m_adrc.rise_time);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Settling Time (s)', m_base.settling_time,m_ieid.settling_time,m_adrc.settling_time);
fprintf('%-22s | %-11.2f%% | %-11.2f%% | %-11.2f%%\n', 'Overshoot',   m_base.overshoot,    m_ieid.overshoot,    m_adrc.overshoot);
fprintf('%-22s | %-12.4e | %-12.4e | %-12.4e\n', 'SS Error',          m_base.sse,          m_ieid.sse,          m_adrc.sse);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Control Energy',    m_base.control_energy,m_ieid.control_energy,m_adrc.control_energy);

%% =========================================================================
%% 5. Plots
%% =========================================================================
fs = 12; lw = 1.8;
C_ref  = [0.3  0.3  0.3];
C_base = [0.8  0.1  0.1];
C_ieid = [0.08 0.40 0.75];
C_adrc = [0.10 0.65 0.30];

% — Figure 1: Tracking
figure('Name','Nonlinear Plant — Output Tracking','Color','w','Position',[100 500 880 360]);
plot(res_base.t, res_base.r, '--', 'Color', C_ref,  'LineWidth', 1.2, 'DisplayName','Reference r'); hold on;
plot(res_base.t, res_base.y, '-',  'Color', C_base, 'LineWidth', lw,  'DisplayName','PID-only (no compensation)');
plot(res_ieid.t, res_ieid.y, '-',  'Color', C_ieid, 'LineWidth', lw,  'DisplayName','IEID');
plot(res_adrc.t, res_adrc.y, '--', 'Color', C_adrc, 'LineWidth', lw,  'DisplayName','ADRC');
hold off; grid on;
xlabel('Time (s)','FontSize',fs); ylabel('Neutron density deviation','FontSize',fs);
title('Nonlinear PWR Plant — Tracking: PID-Only vs IEID vs ADRC','FontSize',fs);
legend('Location','southeast','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 60]);

% — Figure 2: Tracking Error
figure('Name','Nonlinear Plant — Tracking Error','Color','w','Position',[100 80 880 360]);
plot(res_base.t, res_base.r - res_base.y, '-',  'Color', C_base, 'LineWidth', lw, 'DisplayName','PID-only'); hold on;
plot(res_ieid.t, res_ieid.r - res_ieid.y, '-',  'Color', C_ieid, 'LineWidth', lw, 'DisplayName','IEID');
plot(res_adrc.t, res_adrc.r - res_adrc.y, '--', 'Color', C_adrc, 'LineWidth', lw, 'DisplayName','ADRC');
yline(0, ':', 'Color', C_ref, 'LineWidth', 1);
hold off; grid on;
xlabel('Time (s)','FontSize',fs); ylabel('Error  e = r - y','FontSize',fs);
title('Nonlinear PWR Plant — Tracking Error','FontSize',fs);
legend('Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 60]);

% — Figure 3: Disturbance Estimation
figure('Name','Nonlinear Plant — Disturbance','Color','w','Position',[1000 500 880 360]);
plot(res_ieid.t,  res_ieid.d,     '-',       'Color', C_ref,  'LineWidth', 1.2, 'DisplayName','Actual d(t)'); hold on;
plot(res_ieid.t,  res_ieid.d_est, '--',      'Color', C_ieid, 'LineWidth', lw,  'DisplayName','IEID estimate \tilde{d}');
plot(res_adrc.t,  res_adrc.d_est, ':',       'Color', C_adrc, 'LineWidth', lw+0.3, 'DisplayName','ADRC (ESO) z_3/b_0');
hold off; grid on;
xlabel('Time (s)','FontSize',fs); ylabel('Disturbance (equiv. input)','FontSize',fs);
title('Nonlinear PWR Plant — Disturbance Estimation','FontSize',fs);
legend('Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 60]);

% — Figure 4: Control Effort
figure('Name','Nonlinear Plant — Control Effort','Color','w','Position',[1000 80 880 360]);
plot(res_base.t, res_base.u, '-',  'Color', C_base, 'LineWidth', lw, 'DisplayName','PID-only'); hold on;
plot(res_ieid.t, res_ieid.u, '-',  'Color', C_ieid, 'LineWidth', lw, 'DisplayName','IEID');
plot(res_adrc.t, res_adrc.u, '--', 'Color', C_adrc, 'LineWidth', lw, 'DisplayName','ADRC');
hold off; grid on;
xlabel('Time (s)','FontSize',fs); ylabel('Control input u','FontSize',fs);
title('Nonlinear PWR Plant — Control Effort','FontSize',fs);
legend('Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 60]);

fprintf('\nDone. 4 figures generated for nonlinear plant.\n');
