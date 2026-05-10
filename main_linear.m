%% MAIN_LINEAR
% Simulation of the OLD LINEARIZED PWR plant.
% Compares:
%   1. PID-only (with disturbance, no compensation) — baseline
%   2. IEID (disturbance estimation + rejection)

clear; clc; close all;

addpath('plant');
addpath('ieid');
addpath('adrc');
addpath('shared');

%% Parameters
p = get_linear_params();

tspan = [0 60];
opts  = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Initial conditions: deviation form, all zeros (at equilibrium)
x0 = zeros(5, 1);

%% =========================================================================
%% 1. PID-Only Baseline (disturbance present, no compensation)
%% =========================================================================
fprintf('Running PID-Only baseline (old linear plant)...\n');
z0_base = [x0; 0; 0];     % [plant(5); x_int; x_d]

[t_base, z_base] = ode15s(@(t,z) run_pid_only_linear(t, z, p), tspan, z0_base, opts);

res_base.t = t_base;
res_base.y = (p.C * z_base(:, 1:5)')';    % neutron density deviation
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
fprintf('Running IEID (old linear plant)...\n');

% Build param struct fields required by run_closed_loop
p.A_obs  = p.A_real - p.L_obs * p.C;
p.B_obs  = [p.B_real, p.L_obs];
p.L      = p.L_obs;
p.tau    = 0.02;

z0_ieid = zeros(14, 1);   % [plant(5); obs(5); x_int; x_d; lpf(2)]

[t_ieid, z_ieid] = ode15s(@(t,z) run_closed_loop(t, z, p), tspan, z0_ieid, opts);

res_ieid.t     = t_ieid;
res_ieid.y     = (p.C * z_ieid(:, 1:5)')';
res_ieid.r     = 0.1 * (t_ieid >= 20);
res_ieid.d     = arrayfun(@disturbance_fcn, t_ieid);
res_ieid.d_est = z_ieid(:, 13);    % d_tilde (LPF output)
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
%% 3. Metrics
%% =========================================================================
m_base = metrics(t_base, res_base.y, res_base.r, res_base.u);
m_ieid = metrics(t_ieid, res_ieid.y, res_ieid.r, res_ieid.u);

fprintf('\n=== Performance Comparison — Old Linear Plant ===\n');
fprintf('%-22s | %-12s | %-12s\n', 'Metric', 'PID-Only', 'IEID');
fprintf('%s\n', repmat('-', 1, 52));
fprintf('%-22s | %-12.4f | %-12.4f\n', 'RMS Error',         m_base.rms_error,    m_ieid.rms_error);
fprintf('%-22s | %-12.4f | %-12.4f\n', 'Rise Time (s)',      m_base.rise_time,    m_ieid.rise_time);
fprintf('%-22s | %-12.4f | %-12.4f\n', 'Settling Time (s)',  m_base.settling_time,m_ieid.settling_time);
fprintf('%-22s | %-11.2f%% | %-11.2f%%\n', 'Overshoot',      m_base.overshoot,    m_ieid.overshoot);
fprintf('%-22s | %-12.4e | %-12.4e\n', 'SS Error',           m_base.sse,          m_ieid.sse);
fprintf('%-22s | %-12.4f | %-12.4f\n', 'Control Energy',     m_base.control_energy,m_ieid.control_energy);

%% =========================================================================
%% 4. Plots
%% =========================================================================
fs = 12; lw = 1.8;
C_ref  = [0.3  0.3  0.3];
C_base = [0.8  0.1  0.1];
C_ieid = [0.08 0.40 0.75];

% — Figure 1: Tracking
figure('Name','Linear Plant — Output Tracking','Color','w','Position',[100 500 880 360]);
plot(res_base.t, res_base.r, '--', 'Color', C_ref,  'LineWidth', 1.2, 'DisplayName','Reference r'); hold on;
plot(res_base.t, res_base.y, '-',  'Color', C_base, 'LineWidth', lw,  'DisplayName','PID-only (no compensation)');
plot(res_ieid.t, res_ieid.y, '-',  'Color', C_ieid, 'LineWidth', lw,  'DisplayName','IEID');
hold off; grid on;
xlabel('Time (s)','FontSize',fs); ylabel('Neutron density deviation','FontSize',fs);
title('Old Linear Plant — Tracking: PID-Only vs IEID','FontSize',fs);
legend('Location','southeast','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 60]);

% — Figure 2: Tracking Error
figure('Name','Linear Plant — Tracking Error','Color','w','Position',[100 80 880 360]);
plot(res_base.t, res_base.r - res_base.y, '-', 'Color', C_base, 'LineWidth', lw, 'DisplayName','PID-only'); hold on;
plot(res_ieid.t, res_ieid.r - res_ieid.y, '-', 'Color', C_ieid, 'LineWidth', lw, 'DisplayName','IEID');
yline(0, ':', 'Color', C_ref, 'LineWidth', 1);
hold off; grid on;
xlabel('Time (s)','FontSize',fs); ylabel('Error  e = r - y','FontSize',fs);
title('Old Linear Plant — Tracking Error','FontSize',fs);
legend('Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 60]);

% — Figure 3: Disturbance Estimation
figure('Name','Linear Plant — Disturbance','Color','w','Position',[1000 500 880 360]);
plot(res_ieid.t, res_ieid.d,     '-',  'Color', C_ref,  'LineWidth', 1.2, 'DisplayName','Actual d(t)'); hold on;
plot(res_ieid.t, res_ieid.d_est, '--', 'Color', C_ieid, 'LineWidth', lw,  'DisplayName','IEID estimate \tilde{d}');
hold off; grid on;
xlabel('Time (s)','FontSize',fs); ylabel('Disturbance (equiv. input)','FontSize',fs);
title('Old Linear Plant — Disturbance Estimation (IEID)','FontSize',fs);
legend('Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 60]);

% — Figure 4: Control Effort
figure('Name','Linear Plant — Control Effort','Color','w','Position',[1000 80 880 360]);
plot(res_base.t, res_base.u, '-', 'Color', C_base, 'LineWidth', lw, 'DisplayName','PID-only'); hold on;
plot(res_ieid.t, res_ieid.u, '-', 'Color', C_ieid, 'LineWidth', lw, 'DisplayName','IEID');
hold off; grid on;
xlabel('Time (s)','FontSize',fs); ylabel('Control input u','FontSize',fs);
title('Old Linear Plant — Control Effort','FontSize',fs);
legend('Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 60]);

fprintf('\nDone. 4 figures generated for old linear plant.\n');
