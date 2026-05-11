function plotting(res_ieid, res_adrc, res_base_nl, res_ieid_lin, res_base_lin)
%PLOTTING  Generates all comparison plots.
%
%   INPUTS:
%     res_ieid     - IEID results on nonlinear plant
%     res_adrc     - ADRC results on nonlinear plant
%     res_base_nl  - PID-only (no compensation) baseline on nonlinear plant
%     res_ieid_lin - IEID results on old linearized plant
%     res_base_lin - PID-only (no compensation) baseline on old linear plant

fs = 14;
fs_title = 16;
lw = 2.0;

% Colour palette
C_ref   = [0.3 0.3 0.3];      % Dark Grey  — Reference
C_ieid  = [0 0.4470 0.7410];  % Blue       — IEID
C_adrc  = [0.8500 0.3250 0.0980]; % Red    — ADRC

%% =========================================================================
%% Figure 1 — Nonlinear Plant: Tracking with and without compensation
%% =========================================================================
figure('Name','Nonlinear Plant: Tracking Comparison','Color','w', ...
       'Position',[50 550 900 380]);
plot(res_ieid.t,    res_ieid.y,    '-',  'Color', C_ieid, 'LineWidth', lw,  ...
     'DisplayName', 'IEID'); hold on;
plot(res_adrc.t,    res_adrc.y,    '-', 'Color', C_adrc, 'LineWidth', lw,  ...
     'DisplayName', 'ADRC');
plot(res_base_nl.t, res_base_nl.r, '--', 'Color', C_ref,  'LineWidth', 1.5, ...
     'DisplayName', 'Reference r');
hold off;
grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs);
xlabel('Time (s)', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Neutron density deviation', 'FontSize', fs, 'FontWeight', 'bold');
title('Nonlinear Plant — Output Tracking: IEID vs ADRC', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', fs);
xlim([0 60]);

%% =========================================================================
%% Figure 2 — Old Linear Plant: Tracking with and without compensation
%% =========================================================================
figure('Name','Old Linear Plant: Tracking Comparison','Color','w', ...
       'Position',[50 100 900 380]);
plot(res_ieid_lin.t, res_ieid_lin.y, '-',  'Color', C_ieid, 'LineWidth', lw,  ...
     'DisplayName', 'IEID'); hold on;
plot(res_base_lin.t, res_base_lin.r, '--', 'Color', C_ref,  'LineWidth', 1.5, ...
     'DisplayName', 'Reference r');
hold off;
grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs);
xlabel('Time (s)', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Neutron density deviation', 'FontSize', fs, 'FontWeight', 'bold');
title('Old Linear Plant — Output Tracking: IEID', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', fs);
xlim([0 60]);

%% =========================================================================
%% Figure 3 — Nonlinear Plant: Tracking Error
%% =========================================================================
figure('Name','Nonlinear Plant: Tracking Error','Color','w', ...
       'Position',[980 550 900 380]);
plot(res_ieid.t,    res_ieid.r    - res_ieid.y,    '-', 'Color', C_ieid, 'LineWidth', lw, ...
     'DisplayName', 'IEID'); hold on;
plot(res_adrc.t,    res_adrc.r    - res_adrc.y,    '-','Color', C_adrc, 'LineWidth', lw, ...
     'DisplayName', 'ADRC');
yline(0, '--', 'Color', C_ref, 'LineWidth', 1.5);
hold off;
grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs);
xlabel('Time (s)', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Error  e = r - y', 'FontSize', fs, 'FontWeight', 'bold');
title('Nonlinear Plant — Tracking Error', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', fs);
xlim([0 60]);

%% =========================================================================
%% Figure 4 — Old Linear Plant: Tracking Error
%% =========================================================================
figure('Name','Old Linear Plant: Tracking Error','Color','w', ...
       'Position',[980 100 900 380]);
plot(res_ieid_lin.t, res_ieid_lin.r - res_ieid_lin.y, '-', 'Color', C_ieid, 'LineWidth', lw, ...
     'DisplayName', 'IEID'); hold on;
yline(0, '--', 'Color', C_ref, 'LineWidth', 1.5);
hold off;
grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs);
xlabel('Time (s)', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Error  e = r - y', 'FontSize', fs, 'FontWeight', 'bold');
title('Old Linear Plant — Tracking Error', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', fs);
xlim([0 60]);

%% =========================================================================
%% Figure 5 — Disturbance Estimation (Nonlinear Plant)
%% =========================================================================
figure('Name','Disturbance Estimation (Nonlinear)','Color','w', ...
       'Position',[50 300 900 340]);
plot(res_ieid.t,  res_ieid.d_est, '-', 'Color', C_ieid, 'LineWidth', lw,  ...
     'DisplayName', 'IEID estimate \tilde{d}'); hold on;
plot(res_adrc.t,  res_adrc.d_est, '-',  'Color', C_adrc, 'LineWidth', lw,  ...
     'DisplayName', 'ADRC (ESO) z_3/b_0');
plot(res_ieid.t,  res_ieid.d,     '-',  'Color', C_ref,  'LineWidth', 1.5, ...
     'DisplayName', 'Actual d(t)');
hold off;
grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs);
xlabel('Time (s)', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Disturbance (equivalent input)', 'FontSize', fs, 'FontWeight', 'bold');
title('Nonlinear Plant — Disturbance Estimation', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', fs);
xlim([0 60]);

%% =========================================================================
%% Figure 6 — Disturbance Estimation (Old Linear Plant)
%% =========================================================================
figure('Name','Disturbance Estimation (Linear)','Color','w', ...
       'Position',[980 300 900 340]);
plot(res_ieid_lin.t, res_ieid_lin.d_est, '-', 'Color', C_ieid, 'LineWidth', lw,  ...
     'DisplayName', 'IEID estimate \tilde{d}'); hold on;
plot(res_ieid_lin.t, res_ieid_lin.d,     '-',  'Color', C_ref,  'LineWidth', 1.5, ...
     'DisplayName', 'Actual d(t)');
hold off;
grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs);
xlabel('Time (s)', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Disturbance (equivalent input)', 'FontSize', fs, 'FontWeight', 'bold');
title('Old Linear Plant — Disturbance Estimation', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', fs);
xlim([0 60]);

%% =========================================================================
%% Figure 7 — Control Effort Comparison (Nonlinear Plant)
%% =========================================================================
figure('Name','Control Effort (Nonlinear)','Color','w', ...
       'Position',[50 50 900 320]);
plot(res_ieid.t,    res_ieid.u,    '-', 'Color', C_ieid, 'LineWidth', lw, ...
     'DisplayName', 'IEID'); hold on;
plot(res_adrc.t,    res_adrc.u,    '-','Color', C_adrc, 'LineWidth', lw, ...
     'DisplayName', 'ADRC');
hold off;
grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs);
xlabel('Time (s)', 'FontSize', fs, 'FontWeight', 'bold');
ylabel('Control input u (rod reactivity rate)', 'FontSize', fs, 'FontWeight', 'bold');
title('Nonlinear Plant — Control Effort', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location', 'best', 'FontSize', fs);
xlim([0 60]);

end
