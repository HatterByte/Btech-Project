function plotting(res_ieid, res_adrc, res_base_nl, res_ieid_lin, res_base_lin)
%PLOTTING  Generates all comparison plots.
%
%   INPUTS:
%     res_ieid     - IEID results on nonlinear plant
%     res_adrc     - ADRC results on nonlinear plant
%     res_base_nl  - PID-only (no compensation) baseline on nonlinear plant
%     res_ieid_lin - IEID results on old linearized plant
%     res_base_lin - PID-only (no compensation) baseline on old linear plant

fs = 12;
lw = 1.8;

% Colour palette
C_ref   = [0.3  0.3  0.3 ];   % dark grey  — reference
C_base  = [0.8  0.1  0.1 ];   % red        — PID-only (no compensation)
C_ieid  = [0.08 0.40 0.75];   % blue       — IEID
C_adrc  = [0.10 0.65 0.30];   % green      — ADRC

%% =========================================================================
%% Figure 1 — Nonlinear Plant: Tracking with and without compensation
%% =========================================================================
figure('Name','Nonlinear Plant: Tracking Comparison','Color','w', ...
       'Position',[50 550 900 380]);
plot(res_base_nl.t, res_base_nl.r, '--', 'Color', C_ref,  'LineWidth', 1.2, ...
     'DisplayName', 'Reference r'); hold on;
plot(res_base_nl.t, res_base_nl.y, '-',  'Color', C_base, 'LineWidth', lw,  ...
     'DisplayName', 'PID only (no compensation)');
plot(res_ieid.t,    res_ieid.y,    '-',  'Color', C_ieid, 'LineWidth', lw,  ...
     'DisplayName', 'IEID');
plot(res_adrc.t,    res_adrc.y,    '--', 'Color', C_adrc, 'LineWidth', lw,  ...
     'DisplayName', 'ADRC');
hold off;
grid on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('Neutron density deviation', 'FontSize', fs);
title('Nonlinear Plant — Output Tracking: PID-Only vs IEID vs ADRC', 'FontSize', fs);
legend('Location', 'southeast', 'FontSize', fs-1);
set(gca, 'FontSize', fs); xlim([0 60]);

%% =========================================================================
%% Figure 2 — Old Linear Plant: Tracking with and without compensation
%% =========================================================================
figure('Name','Old Linear Plant: Tracking Comparison','Color','w', ...
       'Position',[50 100 900 380]);
plot(res_base_lin.t, res_base_lin.r, '--', 'Color', C_ref,  'LineWidth', 1.2, ...
     'DisplayName', 'Reference r'); hold on;
plot(res_base_lin.t, res_base_lin.y, '-',  'Color', C_base, 'LineWidth', lw,  ...
     'DisplayName', 'PID only (no compensation)');
plot(res_ieid_lin.t, res_ieid_lin.y, '-',  'Color', C_ieid, 'LineWidth', lw,  ...
     'DisplayName', 'IEID');
hold off;
grid on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('Neutron density deviation', 'FontSize', fs);
title('Old Linear Plant — Output Tracking: PID-Only vs IEID', 'FontSize', fs);
legend('Location', 'southeast', 'FontSize', fs-1);
set(gca, 'FontSize', fs); xlim([0 60]);

%% =========================================================================
%% Figure 3 — Nonlinear Plant: Tracking Error
%% =========================================================================
figure('Name','Nonlinear Plant: Tracking Error','Color','w', ...
       'Position',[980 550 900 380]);
plot(res_base_nl.t, res_base_nl.r - res_base_nl.y, '-', 'Color', C_base, 'LineWidth', lw, ...
     'DisplayName', 'PID only'); hold on;
plot(res_ieid.t,    res_ieid.r    - res_ieid.y,    '-', 'Color', C_ieid, 'LineWidth', lw, ...
     'DisplayName', 'IEID');
plot(res_adrc.t,    res_adrc.r    - res_adrc.y,    '--','Color', C_adrc, 'LineWidth', lw, ...
     'DisplayName', 'ADRC');
yline(0, ':', 'Color', C_ref, 'LineWidth', 1);
hold off;
grid on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('Error  e = r - y', 'FontSize', fs);
title('Nonlinear Plant — Tracking Error', 'FontSize', fs);
legend('Location', 'best', 'FontSize', fs-1);
set(gca, 'FontSize', fs); xlim([0 60]);

%% =========================================================================
%% Figure 4 — Old Linear Plant: Tracking Error
%% =========================================================================
figure('Name','Old Linear Plant: Tracking Error','Color','w', ...
       'Position',[980 100 900 380]);
plot(res_base_lin.t, res_base_lin.r - res_base_lin.y, '-', 'Color', C_base, 'LineWidth', lw, ...
     'DisplayName', 'PID only'); hold on;
plot(res_ieid_lin.t, res_ieid_lin.r - res_ieid_lin.y, '-', 'Color', C_ieid, 'LineWidth', lw, ...
     'DisplayName', 'IEID');
yline(0, ':', 'Color', C_ref, 'LineWidth', 1);
hold off;
grid on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('Error  e = r - y', 'FontSize', fs);
title('Old Linear Plant — Tracking Error', 'FontSize', fs);
legend('Location', 'best', 'FontSize', fs-1);
set(gca, 'FontSize', fs); xlim([0 60]);

%% =========================================================================
%% Figure 5 — Disturbance Estimation (Nonlinear Plant)
%% =========================================================================
figure('Name','Disturbance Estimation (Nonlinear)','Color','w', ...
       'Position',[50 300 900 340]);
plot(res_ieid.t,  res_ieid.d,     '-',  'Color', C_ref,  'LineWidth', 1.2, ...
     'DisplayName', 'Actual d(t)'); hold on;
plot(res_ieid.t,  res_ieid.d_est, '--', 'Color', C_ieid, 'LineWidth', lw,  ...
     'DisplayName', 'IEID estimate \tilde{d}');
plot(res_adrc.t,  res_adrc.d_est, ':',  'Color', C_adrc, 'LineWidth', lw+0.3,  ...
     'DisplayName', 'ADRC (ESO) z_3/b_0');
hold off;
grid on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('Disturbance (equivalent input)', 'FontSize', fs);
title('Nonlinear Plant — Disturbance Estimation', 'FontSize', fs);
legend('Location', 'best', 'FontSize', fs-1);
set(gca, 'FontSize', fs); xlim([0 60]);

%% =========================================================================
%% Figure 6 — Disturbance Estimation (Old Linear Plant)
%% =========================================================================
figure('Name','Disturbance Estimation (Linear)','Color','w', ...
       'Position',[980 300 900 340]);
plot(res_ieid_lin.t, res_ieid_lin.d,     '-',  'Color', C_ref,  'LineWidth', 1.2, ...
     'DisplayName', 'Actual d(t)'); hold on;
plot(res_ieid_lin.t, res_ieid_lin.d_est, '--', 'Color', C_ieid, 'LineWidth', lw,  ...
     'DisplayName', 'IEID estimate \tilde{d}');
hold off;
grid on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('Disturbance (equivalent input)', 'FontSize', fs);
title('Old Linear Plant — Disturbance Estimation', 'FontSize', fs);
legend('Location', 'best', 'FontSize', fs-1);
set(gca, 'FontSize', fs); xlim([0 60]);

%% =========================================================================
%% Figure 7 — Control Effort Comparison (Nonlinear Plant)
%% =========================================================================
figure('Name','Control Effort (Nonlinear)','Color','w', ...
       'Position',[50 50 900 320]);
plot(res_base_nl.t, res_base_nl.u, '-', 'Color', C_base, 'LineWidth', lw, ...
     'DisplayName', 'PID only'); hold on;
plot(res_ieid.t,    res_ieid.u,    '-', 'Color', C_ieid, 'LineWidth', lw, ...
     'DisplayName', 'IEID');
plot(res_adrc.t,    res_adrc.u,    '--','Color', C_adrc, 'LineWidth', lw, ...
     'DisplayName', 'ADRC');
hold off;
grid on;
xlabel('Time (s)', 'FontSize', fs);
ylabel('Control input u (rod reactivity rate)', 'FontSize', fs);
title('Nonlinear Plant — Control Effort', 'FontSize', fs);
legend('Location', 'best', 'FontSize', fs-1);
set(gca, 'FontSize', fs); xlim([0 60]);

end
