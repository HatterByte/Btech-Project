%% VALIDATE_ADRC_PAPER
% Validates ADRC implementation against:
%   Ahmad et al. (2023), Annals of Nuclear Energy, 189, 109845.
%
% Produces 3 figures matching paper layout:
%   Fig 1 - Scenario 1 (power tracking): subplot(3,2)
%   Fig 2 - Scenario 2 (Te disturbance): subplot(2,2)
%   Fig 3 - Scenario 3 (rho disturbance): subplot(2,2)
%
% Run: >> validate_adrc_paper

clear; clc; close all;

addpath('plant');
addpath('adrc');
addpath('shared');

%% =========================================================================
%% Parameters
%% =========================================================================
p = get_paper_params();

fprintf('ADRC Validation — Ahmad et al. (2023)  [confirmed from Simulink source]\n');
fprintf('  n0=%.1f | Tf0=%.0f°C | Tl0=%.0f°C\n', p.n0, p.Tf0, p.Tl0);
fprintf('  bo=%-5g  wo=%-5g  wc=%-3g  k1=%-4g  k2=%-4g\n', p.b0, p.wo, p.wc, p.k1, p.k2);

%% =========================================================================
%% Initial Conditions  (100% power equilibrium)
%% =========================================================================
% States: [n; c; Tf; Tl; rho_rod;  z1; z2; z3]
x0_plant = [p.n0; p.n0; p.Tf0; p.Tl0; 0];   % plant at equilibrium
z0_eso   = [p.n0; 0; 0];                       % ESO pre-initialised at y=n0
z0       = [x0_plant; z0_eso];

%% =========================================================================
%% Solve ODE (0–400 s)
%% =========================================================================
tspan = [0 400];
opts  = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'MaxStep', 0.05);

fprintf('Running simulation (0–400 s)...\n');
tic;
[t, Z] = ode15s(@(t,z) run_adrc_paper(t, z, p), tspan, z0, opts);
fprintf('Done in %.1f s  (%d points)\n\n', toc, length(t));

%% =========================================================================
%% Extract signals
%% =========================================================================
n_r    = Z(:,1);
Tf     = Z(:,3);   %#ok
Tl     = Z(:,4);   %#ok
rho_rod= Z(:,5);   %#ok
z1_hat = Z(:,6);
z2_hat = Z(:,7);
z3_hat = Z(:,8);

% Reference signal z_r (paper notation) at each time point
z_r = zeros(size(t));
for k = 1:numel(t)
    tk = t(k);
    if     tk < 50;   z_r(k) = 1.0;
    elseif tk < 100;  z_r(k) = 0.9;
    elseif tk < 150;  z_r(k) = 0.7;
    elseif tk < 156;  z_r(k) = 0.7 + 0.05*(tk - 150);
    else;             z_r(k) = 1.0;
    end
end

% Control input u at each time point
u_sig = zeros(size(t));
for k = 1:numel(t)
    u_bar_k  = p.k1*(z_r(k) - z1_hat(k)) - p.k2*z2_hat(k);
    u_sig(k) = max(p.u_min, min(p.u_max, (u_bar_k - z3_hat(k)) / p.b0));
end

%% =========================================================================
%% Helper: restrict to time window
%% =========================================================================
win = @(ta, tb) (t >= ta & t <= tb);

%% =========================================================================
%% Appearance settings
%% =========================================================================
lw_ref  = 2.0;
lw_out  = 2.0;
fs_ax   = 12;
fs_title= 14;
C_ref   = [0.3 0.3 0.3];    % Dark Grey — reference z_r
C_out   = [0 0.4470 0.7410];% Blue  — plant output n_r
C_ctrl  = [0.8500 0.3250 0.0980]; % Orange-Red — control u

%% =========================================================================
%% FIGURE 1 — Scenario 1: Power Tracking
%%   subplot(3,2): left = n_r, right = u
%%   Row 1: S1a  50–90 s   (step 1.0 → 0.9)
%%   Row 2: S1b  100–140 s (step 0.9 → 0.7)
%%   Row 3: S1c  150–200 s (ramp 0.7 → 1.0)
%% =========================================================================
figure('Name','Fig 1 — Scenario 1: Power Tracking', ...
       'Color','w','Position',[50 50 820 700]);

windows1 = {[50 90], [100 140], [150 200]};
ylims_n1 = {[0.87 1.03], [0.68 0.93], [0.68 1.03]};
labels1  = {'S1a: step 1.0 \rightarrow 0.9', ...
            'S1b: step 0.9 \rightarrow 0.7', ...
            'S1c: ramp \rightarrow 1.0'};

for row = 1:3
    ta = windows1{row}(1);
    tb = windows1{row}(2);
    idx = win(ta, tb);

    % Left: neutron power
    subplot(3,2, 2*row-1);
    plot(t(idx), z_r(idx),  '--', 'Color', C_ref, 'LineWidth', lw_ref, ...
         'DisplayName', 'z_r'); hold on;
    plot(t(idx), n_r(idx),  '-',  'Color', C_out, 'LineWidth', lw_out, ...
         'DisplayName', 'n_r');
    hold off;
    ylim(ylims_n1{row});
    xlim([ta tb]);
    grid on;
    set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
    xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    ylabel('n_r', 'FontSize', fs_ax, 'FontWeight', 'bold');
    title(labels1{row}, 'FontSize', fs_title, 'FontWeight', 'bold');
    if row == 1; legend('z_r','n_r','Location','best','FontSize',fs_ax); end

    % Right: control input u
    subplot(3,2, 2*row);
    plot(t(idx), u_sig(idx), '-', 'Color', C_ctrl, 'LineWidth', lw_out);
    yline(p.u_max, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    yline(p.u_min, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    ylim([p.u_min*1.4  p.u_max*1.4]);
    xlim([ta tb]);
    grid on;
    set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
    xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    ylabel('u  (\deltak/k/s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    title([labels1{row} ' — control'], 'FontSize', fs_title, 'FontWeight', 'bold');
end
sgtitle('Figure 1 — Scenario 1: Power Tracking (compare with paper Fig.)', ...
        'FontSize', 16, 'FontWeight','bold');

%% =========================================================================
%% FIGURE 2 — Scenario 2: T_e (Coolant Inlet Temp) Disturbance
%%   subplot(2,2): left = n_r, right = u
%%   Row 1: S2a  200–230 s  (Te step −2°C)
%%   Row 2: S2b  240–280 s  (Te ramp +7°C)
%% =========================================================================
figure('Name','Fig 2 — Scenario 2: Te Disturbance', ...
       'Color','w','Position',[50 50 820 500]);

windows2 = {[200 230], [240 280]};
ylims_n2 = {[0.97 1.015], [0.97 1.015]};
labels2  = {'S2a: T_e step −2°C', 'S2b: T_e ramp +7°C'};

for row = 1:2
    ta = windows2{row}(1);
    tb = windows2{row}(2);
    idx = win(ta, tb);

    % Left: neutron power
    subplot(2,2, 2*row-1);
    plot(t(idx), z_r(idx), '--', 'Color', C_ref, 'LineWidth', lw_ref, ...
         'DisplayName', 'z_r'); hold on;
    plot(t(idx), n_r(idx), '-',  'Color', C_out, 'LineWidth', lw_out, ...
         'DisplayName', 'n_r');
    hold off;
    ylim(ylims_n2{row});
    xlim([ta tb]);
    grid on;
    set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
    xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    ylabel('n_r', 'FontSize', fs_ax, 'FontWeight', 'bold');
    title(labels2{row}, 'FontSize', fs_title, 'FontWeight', 'bold');
    if row == 1; legend('z_r','n_r','Location','best','FontSize',fs_ax); end

    % Right: control input u
    subplot(2,2, 2*row);
    plot(t(idx), u_sig(idx), '-', 'Color', C_ctrl, 'LineWidth', lw_out);
    yline(p.u_max, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    yline(p.u_min, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    ylim([p.u_min*1.4  p.u_max*1.4]);
    xlim([ta tb]);
    grid on;
    set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
    xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    ylabel('u  (\deltak/k/s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    title([labels2{row} ' — control'], 'FontSize', fs_title, 'FontWeight', 'bold');
end
sgtitle('Figure 2 — Scenario 2: T_e Disturbance (compare with paper Fig.)', ...
        'FontSize', 16, 'FontWeight','bold');

%% =========================================================================
%% FIGURE 3 — Scenario 3: Reactivity Disturbance
%%   subplot(2,2): left = n_r, right = u
%%   Row 1: S3a  290–315 s  (rho step)
%%   Row 2: S3b  320–370 s  (rho ramp)
%% =========================================================================
figure('Name','Fig 3 — Scenario 3: Reactivity Disturbance', ...
       'Color','w','Position',[50 50 820 500]);

windows3 = {[290 315], [320 370]};
ylims_n3 = {[0.97 1.015], [0.97 1.015]};
labels3  = {'S3a: \delta\rho step', 'S3b: \delta\rho ramp'};

for row = 1:2
    ta = windows3{row}(1);
    tb = windows3{row}(2);
    idx = win(ta, tb);

    % Left: neutron power
    subplot(2,2, 2*row-1);
    plot(t(idx), z_r(idx), '--', 'Color', C_ref, 'LineWidth', lw_ref, ...
         'DisplayName', 'z_r'); hold on;
    plot(t(idx), n_r(idx), '-',  'Color', C_out, 'LineWidth', lw_out, ...
         'DisplayName', 'n_r');
    hold off;
    ylim(ylims_n3{row});
    xlim([ta tb]);
    grid on;
    set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
    xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    ylabel('n_r', 'FontSize', fs_ax, 'FontWeight', 'bold');
    title(labels3{row}, 'FontSize', fs_title, 'FontWeight', 'bold');
    if row == 1; legend('z_r','n_r','Location','best','FontSize',fs_ax); end

    % Right: control input u
    subplot(2,2, 2*row);
    plot(t(idx), u_sig(idx), '-', 'Color', C_ctrl, 'LineWidth', lw_out);
    yline(p.u_max, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    yline(p.u_min, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
    ylim([p.u_min*1.4  p.u_max*1.4]);
    xlim([ta tb]);
    grid on;
    set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
    xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    ylabel('u  (\deltak/k/s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
    title([labels3{row} ' — control'], 'FontSize', fs_title, 'FontWeight', 'bold');
end
sgtitle('Figure 3 — Scenario 3: Reactivity Disturbance (compare with paper Fig.)', ...
        'FontSize', 16, 'FontWeight','bold');

%% =========================================================================
%% Steady-state summary
%% =========================================================================
ss = t >= 380;
fprintf('=== Steady-state check (t=380–400 s) ===\n');
fprintf('  n_r  mean = %.6f  (target: 1.0)\n', mean(n_r(ss)));
fprintf('  |n_r - z_r| mean = %.2e\n', mean(abs(n_r(ss) - z_r(ss))));

%% =========================================================================
%% FIGURE 4 — ESO Disturbance Estimation (Suggestion 4)
%% Compare z3_hat against the known "actual disturbance" signal in each scenario
%% =========================================================================

% Reconstruct actual disturbances at each time point
d_rho_sig = zeros(size(t));   % Scenario 3 reactivity disturbance
Te_sig    = p.Te0 * ones(size(t));
for k = 1:numel(t)
    tk = t(k);
    if tk >= 200; Te_sig(k) = Te_sig(k) - 2.0; end
    if tk >= 240; Te_sig(k) = Te_sig(k) + min(2.0*(tk-240), 7.0); end
    if tk >= 290; d_rho_sig(k) = -0.0002; end
    if tk >= 320; d_rho_sig(k) = d_rho_sig(k) + min(0.0002*(tk-320), 0.002); end
end

% The ESO z3 estimates the "total disturbance" in the ẋ2 = z3 + b0*u equation.
% It captures EVERYTHING unexplained: temperature coupling, Te variation, d_rho.
% We normalise z3 back to reactivity units: z3/b0 ≈ equivalent input disturbance.
d_est_norm = z3_hat / p.b0;

figure('Name','Fig 4 — ESO Disturbance Estimation','Color','w', ...
       'Position',[50 50 1100 480]);

subplot(3,1,1);
plot(t, d_rho_sig, '-',  'Color', C_ref, 'LineWidth', 1.5, 'DisplayName','d\rho (actual reactivity disturb)'); hold on;
plot(t, d_est_norm,'-', 'Color', C_out, 'LineWidth', lw_out, 'DisplayName','z_3/b_0 (ESO estimate)');
hold off; grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold'); 
ylabel('δk/k', 'FontSize', fs_ax, 'FontWeight', 'bold');
title('ESO Disturbance Estimate vs Actual  (z_3/b_0 captures total disturbance)', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location','best', 'FontSize', fs_ax); xlim([0 400]);

subplot(3,1,2);
yyaxis left;
plot(t, Te_sig - p.Te0, '-', 'Color', C_ctrl, 'LineWidth', lw_out, 'DisplayName','\DeltaTe (°C)');
ylabel('\DeltaTe from nominal (°C)', 'FontSize', fs_ax, 'FontWeight', 'bold');
yyaxis right;
plot(t, d_est_norm, '-', 'Color', C_out, 'LineWidth', lw_out, 'DisplayName','z_3/b_0');
ylabel('ESO disturbance estimate', 'FontSize', fs_ax, 'FontWeight', 'bold');
grid on; 
set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold');
title('Scenario 2 context: T_e variation vs ESO response', 'FontSize', fs_title, 'FontWeight', 'bold');
xlim([180 300]); 

subplot(3,1,3);
plot(t, d_rho_sig, '-',  'Color', C_ref, 'LineWidth', 1.5, 'DisplayName','d\rho actual'); hold on;
plot(t, d_est_norm,'-', 'Color', C_out, 'LineWidth', lw_out, 'DisplayName','z_3/b_0 estimate');
hold off; grid on;
set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold'); 
ylabel('δk/k', 'FontSize', fs_ax, 'FontWeight', 'bold');
title('Scenario 3 zoom: Reactivity disturbance estimation', 'FontSize', fs_title, 'FontWeight', 'bold');
legend('Location','best', 'FontSize', fs_ax); xlim([280 400]);

sgtitle('Figure 4 — ESO Disturbance Estimation (z_3/b_0) — Key ADRC vs IEID comparison metric', ...
        'FontSize', 16, 'FontWeight','bold');

fprintf('\n4 figures generated.\n');
fprintf('Note: z3/b0 (ESO estimate) vs actual disturbances — compare with IEID d_tilde later.\n');

