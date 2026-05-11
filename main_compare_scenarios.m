%% MAIN_COMPARE_SCENARIOS
% Comprehensive comparison of PID, IEID, ADRC, and CESO across the 
% original ADRC paper scenarios, adapted to the 50% power operating point.

clear; clc; close all;

% Add folders to path
addpath('plant');
addpath('ieid');
addpath('adrc');
addpath('shared');

% ODE options
opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 0.1);

% Base nonlinear parameters
p_base = get_nonlinear_params();

%% Scenarios Definition
scenarios = {
    struct('name', 'Scenario 1: Power Tracking (Steps & Ramps)', ...
           'tspan', [0 100], ...
           'r_dev_fcn', @(t) tracking_ref(t), ...
           'Te_fcn', @(t) 290, ...
           'd_rho_fcn', @(t) 0, ...
           'd_fcn', @(t) disturbance_fcn(t)), ...
           
    struct('name', 'Scenario 2: Coolant Inlet Temp (T_e) Disturbance', ...
           'tspan', [0 90], ...
           'r_dev_fcn', @(t) 0, ...
           'Te_fcn', @(t) te_dist(t), ...
           'd_rho_fcn', @(t) 0, ...
           'd_fcn', @(t) disturbance_fcn(t)), ...
           
    struct('name', 'Scenario 3: External Reactivity (\delta\rho) Disturbance', ...
           'tspan', [0 90], ...
           'r_dev_fcn', @(t) 0, ...
           'Te_fcn', @(t) 290, ...
           'd_rho_fcn', @(t) rho_dist(t), ...
           'd_fcn', @(t) disturbance_fcn(t))
};

%% Run Simulations
results = cell(length(scenarios), 3); % rows: scenarios, cols: [IEID, ESO-ADRC, CESO-ADRC]
controllers = {'IEID', 'ESO-ADRC', 'CESO-ADRC'};

x0_nl = [0.5; 0.5; 650; 314; 0]; % Equilibrium at 50% power

for s = 1:length(scenarios)
    fprintf('\n======================================================\n');
    fprintf('Running %s\n', scenarios{s}.name);
    fprintf('======================================================\n');
    
    p = p_base;
    p.r_dev_fcn = scenarios{s}.r_dev_fcn;
    p.Te_fcn = scenarios{s}.Te_fcn;
    p.d_rho_fcn = scenarios{s}.d_rho_fcn;
    p.d_fcn = scenarios{s}.d_fcn;
    tspan = scenarios{s}.tspan;

    % 1. IEID
    fprintf('  -> IEID...\n');
    z0_ieid = [x0_nl; zeros(5,1); 0; 0; 0; 0];
    [t_ieid, z_ieid] = ode15s(@(t,z) run_ieid_closed_loop(t, z, p), tspan, z0_ieid, opts);
    res_ieid.t = t_ieid;
    res_ieid.y = z_ieid(:,1) - 0.5;
    res_ieid.r = arrayfun(p.r_dev_fcn, t_ieid);
    res_ieid.u = compute_ieid_u(t_ieid, z_ieid, p);
    results{s, 1} = res_ieid;
    
    % 2. ESO-ADRC
    fprintf('  -> ESO-ADRC...\n');
    z0_adrc = [x0_nl; zeros(3,1)];
    [t_adrc, z_adrc] = ode15s(@(t,z) run_adrc_closed_loop(t, z, p), tspan, z0_adrc, opts);
    res_adrc.t = t_adrc;
    res_adrc.y = z_adrc(:,1) - 0.5;
    res_adrc.r = arrayfun(p.r_dev_fcn, t_adrc);
    res_adrc.u = compute_adrc_u(t_adrc, z_adrc, p);
    results{s, 2} = res_adrc;
    
    % 3. CESO-ADRC
    fprintf('  -> CESO-ADRC...\n');
    z0_ceso = [x0_nl; zeros(6,1)];
    [t_ceso, z_ceso] = ode15s(@(t,z) run_ceso_closed_loop(t, z, p), tspan, z0_ceso, opts);
    res_ceso.t = t_ceso;
    res_ceso.y = z_ceso(:,1) - 0.5;
    res_ceso.r = arrayfun(p.r_dev_fcn, t_ceso);
    res_ceso.u = compute_ceso_u(t_ceso, z_ceso, p);
    results{s, 3} = res_ceso;
end

%% Plotting
% Vibrant, distinct colors for thesis
colors = {[0.3 0.3 0.3], ...          % Reference (Dark Grey)
          [0 0.4470 0.7410], ...      % IEID (Blue)
          [0.8500 0.3250 0.0980], ... % ESO-ADRC (Red)
          [0.4660 0.6740 0.1880]};    % CESO-ADRC (Green)
          
lw = 2.0;
fs_title = 14;
fs_ax = 12;

for s = 1:length(scenarios)
    fig = figure('Name', scenarios{s}.name, 'Color', 'w', 'Position', [100+s*50, 100, 1000, 400]);
    
    % Output Tracking only
    for c = [2, 3, 1]
        plot(results{s,c}.t, results{s,c}.y + 0.5, '-', 'Color', colors{c+1}, 'LineWidth', lw, 'DisplayName', controllers{c}); hold on;
    end
    plot(results{s,1}.t, results{s,1}.r + 0.5, '--', 'Color', colors{1}, 'LineWidth', 1.7, 'DisplayName', 'Reference');
    
    hold off; grid on;
    set(gca, 'GridAlpha', 0.15, 'FontSize', fs_ax);
    xlabel('Time (s)', 'FontSize', fs_ax, 'FontWeight', 'bold'); 
    ylabel('Relative Power n_r', 'FontSize', fs_ax, 'FontWeight', 'bold');
    title(scenarios{s}.name, 'FontSize', fs_title, 'FontWeight', 'bold');
    legend('Location', 'best', 'FontSize', fs_ax);
end
fprintf('\nSimulation and plotting complete!\n');

%% Helper Functions for Scenario Signals
function r = tracking_ref(t)
    if t < 25;       r = 0;
    elseif t < 50;   r = 0.05;
    elseif t < 75;   r = -0.05;
    elseif t < 85;   r = -0.05 + 0.01*(t-75);
    else;            r = 0.05;
    end
end

function Te = te_dist(t)
    Te = 290;
    if t >= 20; Te = Te - 2.0; end
    if t >= 60; Te = Te + min(2.0*(t - 60), 7.0); end
end

function d_rho = rho_dist(t)
    d_rho = 0;
    if t >= 20; d_rho = -0.0002; end
    if t >= 60; d_rho = d_rho + min(0.0002*(t - 60), 0.002); end
end

%% Helper Functions for Control Signal Extraction
function u = compute_pid_u(t_vec, z_mat, p)
    u = zeros(length(t_vec), 1);
    for k = 1:length(t_vec)
        y  = z_mat(k,1) - 0.5;
        e  = p.r_dev_fcn(t_vec(k)) - y;
        xi = z_mat(k,6);
        xd = z_mat(k,7);
        u(k) = p.Kp*e + p.Ki*xi + p.Kd*p.N*(e - xd);
    end
end

function u = compute_ieid_u(t_vec, z_mat, p)
    u = zeros(length(t_vec), 1);
    for k = 1:length(t_vec)
        y       = z_mat(k,1) - 0.5;
        e       = p.r_dev_fcn(t_vec(k)) - y;
        xi      = z_mat(k,11);
        xd      = z_mat(k,12);
        dtilde  = z_mat(k,13);
        u_f     = p.Kp*e + p.Ki*xi + p.Kd*p.N*(e - xd);
        u(k)    = u_f - dtilde;
    end
end

function u = compute_adrc_u(t_vec, z_mat, p)
    u = zeros(length(t_vec), 1);
    for k = 1:length(t_vec)
        z_eso = z_mat(k,6:8)';
        r_dev = p.r_dev_fcn(t_vec(k));
        u_bar = p.k1*(r_dev - z_eso(1)) - p.k2*z_eso(2);
        u(k)  = (u_bar - z_eso(3)) / p.b0;
    end
end

function u = compute_ceso_u(t_vec, z_mat, p)
    u = zeros(length(t_vec), 1);
    for k = 1:length(t_vec)
        z_ceso = z_mat(k,6:11)';
        z1 = z_ceso(1:3);
        z2 = z_ceso(4:6);
        x_hat = [z2(1); z2(2); z2(3) + z1(3)];
        r_dev = p.r_dev_fcn(t_vec(k));
        u(k) = adrc_control_law(r_dev, x_hat, p);
    end
end
