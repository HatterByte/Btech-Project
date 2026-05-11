%% MAIN_NONLINEAR
% Simulation of the nonlinear PWR plant.
% Compares:
%   1. IEID disturbance compensation
%   2. ADRC with standard ESO
%   3. ADRC with two-level CESO

clear; clc; close all;

addpath('plant');
addpath('ieid');
addpath('adrc');
addpath('shared');

%% Parameters
p = get_nonlinear_params();
p.r_dev_fcn = @(t) 0.1 * (t >= 25);

tspan = [0 40];
opts  = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% Common initial conditions: [n; c; Tf; Tl; rho_rod]
x0 = [0.5; 0.5; 650; 314; 0];

%% 1. IEID
fprintf('Running IEID (nonlinear plant)...\n');
z0_ieid = [x0; zeros(5,1); 0; 0; 0; 0];

[t_ieid, z_ieid] = ode15s(@(t,z) run_ieid_closed_loop(t, z, p), tspan, z0_ieid, opts);

res_ieid.t     = t_ieid;
res_ieid.y     = z_ieid(:, 1) - 0.5;
res_ieid.r     = 0.1 * (t_ieid >= 25);
res_ieid.d     = arrayfun(@disturbance_fcn, t_ieid);
res_ieid.d_est = z_ieid(:, 13);
res_ieid.u     = zeros(length(t_ieid), 1);
for k = 1:length(t_ieid)
    e      = res_ieid.r(k) - res_ieid.y(k);
    xi     = z_ieid(k, 11);
    xd     = z_ieid(k, 12);
    dtilde = z_ieid(k, 13);
    u_f    = p.Kp*e + p.Ki*xi + p.Kd*p.N*(e - xd);
    res_ieid.u(k) = u_f - dtilde;
end

%% 2. ADRC with standard ESO
fprintf('Running ADRC-ESO (nonlinear plant)...\n');
z0_adrc = [x0; zeros(3,1)];

[t_adrc, z_adrc] = ode15s(@(t,z) run_adrc_closed_loop(t, z, p), tspan, z0_adrc, opts);

res_adrc.t     = t_adrc;
res_adrc.y     = z_adrc(:, 1) - 0.5;
res_adrc.r     = 0.1 * (t_adrc >= 25);
res_adrc.d     = arrayfun(@disturbance_fcn, t_adrc);
res_adrc.d_est = z_adrc(:, 8) / p.b0;
res_adrc.u     = zeros(length(t_adrc), 1);
for k = 1:length(t_adrc)
    r_dev = 0.1 * (t_adrc(k) >= 25);
    x_hat = z_adrc(k, 6:8)';
    res_adrc.u(k) = adrc_control_law(r_dev, x_hat, p);
end

%% 3. ADRC with two-level CESO
fprintf('Running CESO-ADRC (nonlinear plant)...\n');
z0_ceso = [x0; zeros(6,1)];

[t_ceso, z_ceso] = ode15s(@(t,z) run_ceso_closed_loop(t, z, p), tspan, z0_ceso, opts);

res_ceso.t     = t_ceso;
res_ceso.y     = z_ceso(:, 1) - 0.5;
res_ceso.r     = 0.1 * (t_ceso >= 25);
res_ceso.d     = arrayfun(@disturbance_fcn, t_ceso);
res_ceso.d_est = zeros(length(t_ceso), 1);
res_ceso.u     = zeros(length(t_ceso), 1);
for k = 1:length(t_ceso)
    z1 = z_ceso(k, 6:8)';
    z2 = z_ceso(k, 9:11)';
    x_hat = [z2(1); z2(2); z2(3) + z1(3)];
    r_dev = 0.1 * (t_ceso(k) >= 25);
    res_ceso.u(k) = adrc_control_law(r_dev, x_hat, p);
    res_ceso.d_est(k) = x_hat(3) / p.b0;
end

%% Metrics
m_ieid = metrics(t_ieid, res_ieid.y, res_ieid.r, res_ieid.u);
m_adrc = metrics(t_adrc, res_adrc.y, res_adrc.r, res_adrc.u);
m_ceso = metrics(t_ceso, res_ceso.y, res_ceso.r, res_ceso.u);

fprintf('\n=== Performance Comparison: Nonlinear PWR Plant ===\n');
fprintf('%-22s | %-12s | %-12s | %-12s\n', 'Metric', 'IEID', 'ADRC-ESO', 'CESO-ADRC');
fprintf('%s\n', repmat('-', 1, 65));
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'RMS Error',        m_ieid.rms_error,     m_adrc.rms_error,     m_ceso.rms_error);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Rise Time (s)',     m_ieid.rise_time,     m_adrc.rise_time,     m_ceso.rise_time);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Settling Time (s)', m_ieid.settling_time, m_adrc.settling_time, m_ceso.settling_time);
fprintf('%-22s | %-11.2f%% | %-11.2f%% | %-11.2f%%\n', 'Overshoot',   m_ieid.overshoot,     m_adrc.overshoot,     m_ceso.overshoot);
fprintf('%-22s | %-12.4e | %-12.4e | %-12.4e\n', 'SS Error',          m_ieid.sse,           m_adrc.sse,           m_ceso.sse);
fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Control Energy',    m_ieid.control_energy, m_adrc.control_energy, m_ceso.control_energy);

%% Plots
fs = 12; lw = 1.8;
C_ref  = [0.30 0.30 0.30];
C_base = [0.80 0.10 0.10];
C_ieid = [0 0.4470 0.7410];
C_adrc = [0.8500 0.3250 0.0980]; % Red
C_ceso = [0.4660 0.6740 0.1880]; % Green

figure('Name','Nonlinear Plant - Output Tracking','Color','w','Position',[100 500 920 380]);
plot(res_ieid.t, res_ieid.y, '-',  'Color', C_ieid, 'LineWidth', lw,  'DisplayName','IEID'); hold on;
plot(res_adrc.t, res_adrc.y, '-', 'Color', C_adrc, 'LineWidth', lw,  'DisplayName','ADRC-ESO');
plot(res_ceso.t, res_ceso.y, '-', 'Color', C_ceso, 'LineWidth', lw,  'DisplayName','CESO-ADRC');
plot(res_ieid.t, res_ieid.r, '--', 'Color', C_ref,  'LineWidth', 1.7, 'DisplayName','Reference');
hold off; grid on;
xlabel('Time (s)','FontSize',fs);
ylabel('Neutron density deviation','FontSize',fs);
title('Nonlinear PWR Plant - Output Tracking','FontSize',fs);
legend('Location','southeast','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 40]);

figure('Name','Nonlinear Plant - Disturbance Estimation','Color','w','Position',[1000 500 920 380]);
plot(res_ieid.t, res_ieid.d_est, '-', 'Color', C_ieid, 'LineWidth', lw, 'DisplayName','IEID estimate'); hold on;
plot(res_adrc.t, res_adrc.d_est, '-', 'Color', C_adrc, 'LineWidth', lw, 'DisplayName','ADRC-ESO estimate');
plot(res_ceso.t, res_ceso.d_est, '-', 'Color', C_ceso, 'LineWidth', lw, 'DisplayName','CESO-ADRC estimate');
plot(res_ieid.t, res_ieid.d, '-', 'Color', C_ref, 'LineWidth', 1.7, 'DisplayName','Actual d(t)');
hold off; grid on;
xlabel('Time (s)','FontSize',fs);
ylabel('Disturbance (equiv. input)','FontSize',fs);
title('Nonlinear PWR Plant - Disturbance Estimation','FontSize',fs);
legend('Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 40]);

figure('Name','Nonlinear Plant - Control Effort','Color','w','Position',[1000 80 920 380]);
plot(res_ieid.t, res_ieid.u, '-',  'Color', C_ieid, 'LineWidth', lw, 'DisplayName','IEID'); hold on;
plot(res_adrc.t, res_adrc.u, '-', 'Color', C_adrc, 'LineWidth', lw, 'DisplayName','ADRC-ESO');
plot(res_ceso.t, res_ceso.u, '-', 'Color', C_ceso, 'LineWidth', lw, 'DisplayName','CESO-ADRC');
hold off; grid on;
xlabel('Time (s)','FontSize',fs);
ylabel('Control input u','FontSize',fs);
title('Nonlinear PWR Plant - Control Effort','FontSize',fs);
legend('Location','best','FontSize',fs-1);
set(gca,'FontSize',fs); xlim([0 40]);

fprintf('\nDone. 3 figures generated for nonlinear plant.\n');
