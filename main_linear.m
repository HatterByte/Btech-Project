%% MAIN_LINEAR
% Simulation of the OLD LINEARIZED PWR plant with multiple cases:
%   Case 1: Nominal Tracking (No disturbance)
%   Case 2: Disturbance Rejection
%   Case 3: Parameter Uncertainty
%
% Compares:
%   1. PID-only (baseline)
%   2. IEID (disturbance estimation + rejection)
%   3. MPC (Predictive control)

clear; clc; close all;

addpath('plant');
addpath('ieid');
addpath('adrc');
addpath('shared');
addpath('mpc');

%% Parameters
p = get_linear_params();
tspan = [0 60];
opts  = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);
x0    = zeros(5, 1);

% Prepare MPC
mpcobj = setup_mpc(p);

%% =========================================================================
%% RUN CASES
%% =========================================================================

cases = {'Nominal', 'Disturbance', 'Uncertainty'};

for i = 1:length(cases)
    case_name = cases{i};
    fprintf('\n--- RUNNING CASE: %s ---\n', case_name);
    
    % Setup case specific parameters
    p_case = p;
    if strcmp(case_name, 'Nominal')
        % No disturbance function should be used or it should return 0
        dist_func = @(t) 0;
    elseif strcmp(case_name, 'Disturbance')
        dist_func = @disturbance_fcn;
    elseif strcmp(case_name, 'Uncertainty')
        dist_func = @disturbance_fcn;
        % Perturb A matrix by 10%
        p_case.A_real = p.A_real * 1.1; 
        p_case.B_real = p.B_real * 0.9;
    end
    
    % 1. PID-Only
    z0_base = [x0; 0; 0];
    [t_base, z_base] = ode15s(@(t,z) run_pid_linear_wrapper(t, z, p_case, dist_func), tspan, z0_base, opts);
    res(i).base = process_pid_res(t_base, z_base, p_case);
    
    % 2. IEID
    p_ieid = p_case;
    p_ieid.A_obs = p.A_real - p.L_obs * p.C; % Observer uses NOMINAL model
    p_ieid.B_obs = [p.B_real, p.L_obs];      % Observer uses NOMINAL model
    z0_ieid = zeros(14, 1);
    [t_ieid, z_ieid] = ode15s(@(t,z) run_ieid_linear_wrapper(t, z, p_ieid, dist_func), tspan, z0_ieid, opts);
    res(i).ieid = process_ieid_res(t_ieid, z_ieid, p_ieid, dist_func);
    
    % 3. MPC
    [t_mpc, y_mpc, u_mpc] = run_mpc_linear_sim_wrapper(p_case, mpcobj, tspan, x0, dist_func);
    res(i).mpc.t = t_mpc;
    res(i).mpc.y = y_mpc;
    res(i).mpc.u = u_mpc;
    res(i).mpc.r = 0.1 * (t_mpc >= 20);
    
    % Compute Metrics
    res(i).m_base = metrics(t_base, res(i).base.y, res(i).base.r, res(i).base.u);
    res(i).m_ieid = metrics(t_ieid, res(i).ieid.y, res(i).ieid.r, res(i).ieid.u);
    res(i).m_mpc  = metrics(t_mpc,  res(i).mpc.y,  res(i).mpc.r,  res(i).mpc.u);
    
    % Print Results Table
    print_table(case_name, res(i).m_base, res(i).m_ieid, res(i).m_mpc);
end

%% =========================================================================
%% PLOTTING (Case 2: Disturbance Rejection)
%% =========================================================================
plot_comparison(res(2), 'Case 2: Disturbance Rejection');

%% =========================================================================
%% HELPER FUNCTIONS
%% =========================================================================

function dz = run_pid_linear_wrapper(t, z, p, dist_func)
    x_plant = z(1:5); e = (0.1*(t>=20)) - (p.C*x_plant);
    u = p.Kp*e + p.Ki*z(6) + p.Kd*p.N*(e - z(7));
    dx_plant = p.A_real*x_plant + p.B_real*u + p.B_d*dist_func(t);
    dz = [dx_plant; e; -p.N*z(7) + p.N*e];
end

function dz = run_ieid_linear_wrapper(t, z, p, dist_func)
    % Plant states
    x = z(1:5); y = p.C*x; r = 0.1*(t>=20); e = r - y;
    % Observer and Control
    x_hat = z(6:10); x_int = z(11); x_d = z(12); xf1 = z(13); xf2 = z(14);
    u_f = p.Kp*e + p.Ki*x_int + p.Kd*p.N*(e - x_d);
    u = u_f - xf1;
    
    % Dynamics
    dx = p.A_real*x + p.B_real*u + p.B_d*dist_func(t);
    
    % IEID Estimator (standard form)
    y_hat = p.C * x_hat;
    residual = y - y_hat;
    d_hat = pinv(p.B_real) * (p.L_obs * residual) + (u_f - u);
    
    % Observer and Filter
    dx_hat = p.A_obs*x_hat + p.B_obs*[u_f; y];
    tau = 0.02; dxf1 = xf2; dxf2 = (d_hat - xf1 - 1.41*tau*xf2)/(tau^2);
    dz = [dx; dx_hat; e; -p.N*x_d + p.N*e; dxf1; dxf2];
end

function [t_out, y_out, u_out] = run_mpc_linear_sim_wrapper(p, mpcobj, tspan, x0, dist_func)
    Ts = mpcobj.Ts; t_steps = tspan(1):Ts:tspan(2);
    x = x0; t_out = []; y_out = []; u_out = []; 
    x_mpc_state = mpcstate(mpcobj); % Renamed variable to avoid shadowing function
    for k = 1:length(t_steps)-1
        y = p.C*x; r = 0.1*(t_steps(k)>=20);
        u = mpcmove(mpcobj, x_mpc_state, y, r);
        [t_ode, x_ode] = ode15s(@(t,x) p.A_real*x + p.B_real*u + p.B_d*dist_func(t), [t_steps(k) t_steps(k+1)], x);
        t_out = [t_out; t_ode]; y_out = [y_out; (p.C*x_ode')']; u_out = [u_out; u*ones(length(t_ode),1)];
        x = x_ode(end,:)';
    end
end

function res = process_pid_res(t, z, p)
    res.t = t; res.y = (p.C*z(:,1:5)')'; res.r = 0.1*(t>=20);
    e = res.r - res.y; res.u = p.Kp*e + p.Ki*z(:,6) + p.Kd*p.N*(e - z(:,7));
end

function res = process_ieid_res(t, z, p, dist_func)
    res.t = t; res.y = (p.C*z(:,1:5)')'; res.r = 0.1*(t>=20);
    e = res.r - res.y; u_f = p.Kp*e + p.Ki*z(:,11) + p.Kd*p.N*(e - z(:,12));
    res.u = u_f - z(:,13); res.d = arrayfun(dist_func, t); res.d_est = z(:,13);
end

function print_table(name, m1, m2, m3)
    fprintf('\n=== %s Performance ===\n', name);
    fprintf('%-22s | %-12s | %-12s | %-12s\n', 'Metric', 'PID-Only', 'IEID', 'MPC');
    fprintf('%s\n', repmat('-', 1, 65));
    fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'RMS Error', m1.rms_error, m2.rms_error, m3.rms_error);
    fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Rise Time (s)', m1.rise_time, m2.rise_time, m3.rise_time);
    fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Overshoot (%)', m1.overshoot, m2.overshoot, m3.overshoot);
    fprintf('%-22s | %-12.4f | %-12.4f | %-12.4f\n', 'Control Energy', m1.control_energy, m2.control_energy, m3.control_energy);
end

function plot_comparison(res, title_str)
    fs = 12; lw = 1.8; C_ref = [0.3 0.3 0.3]; C_base = [0.8 0.1 0.1]; C_ieid = [0.08 0.40 0.75]; C_mpc = [0.10 0.65 0.30];
    figure('Name',title_str,'Color','w','Position', [100 100 1000 700]);
    subplot(2,1,1);
    plot(res.base.t, res.base.r, '--', 'Color', C_ref, 'LineWidth', 1.2, 'DisplayName', 'Reference'); hold on;
    plot(res.base.t, res.base.y, '-', 'Color', C_base, 'LineWidth', lw, 'DisplayName', 'PID-Only');
    plot(res.ieid.t, res.ieid.y, '-', 'Color', C_ieid, 'LineWidth', lw, 'DisplayName', 'IEID');
    plot(res.mpc.t, res.mpc.y, '--', 'Color', C_mpc, 'LineWidth', lw, 'DisplayName', 'MPC');
    grid on; ylabel('Output Deviation'); title(title_str); legend('Location','best');
    subplot(2,1,2);
    plot(res.base.t, res.base.u, '-', 'Color', C_base, 'LineWidth', lw, 'DisplayName', 'PID-Only'); hold on;
    plot(res.ieid.t, res.ieid.u, '-', 'Color', C_ieid, 'LineWidth', lw, 'DisplayName', 'IEID');
    plot(res.mpc.t, res.mpc.u, '--', 'Color', C_mpc, 'LineWidth', lw, 'DisplayName', 'MPC');
    grid on; ylabel('Control Effort (u)'); xlabel('Time (s)');
end
