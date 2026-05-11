%% RUN_IEID_OPT_NL
% Optimizes IEID parameters for the NONLINEAR PWR plant using GCRA.
% Scaled down to 30 agents and 20 iterations due to nonlinear ODE stiffness.

clear; clc;


addpath('opt');
addpath('opt/core');
addpath('opt/obj_function');
addpath('shared');
addpath('plant');
addpath('ieid');

%% ---- GCRA Settings ----
SearchAgents_no = 30;
Max_iteration   = 20;

%% ---- Search bounds ----
%        Kp     Ki      Kd      tau     log10(Q1)
lb = [  1.0,   0.1,   0.001,  0.003,      -1  ];
ub = [ 25.0,  20.0,   0.200,  0.100,       6  ];
dim = 5;

% The objective function handle
fobj = @ieid_objective_nl;

fprintf('========================================================\n');
fprintf('=== IEID NONLINEAR Parameter Optimization (GCRA)     ===\n');
fprintf('========================================================\n');
fprintf('Agents: %d | Iterations: %d | Dimensions: %d\n', SearchAgents_no, Max_iteration, dim);
fprintf('Bounds:\n');
fprintf('  Kp  : [%.1f , %.1f]\n', lb(1), ub(1));
fprintf('  Ki  : [%.1f , %.1f]\n', lb(2), ub(2));
fprintf('  Kd  : [%.3f , %.3f]\n', lb(3), ub(3));
fprintf('  tau : [%.4f, %.3f] s\n', lb(4), ub(4));
fprintf('  Q1  : [1e%d, 1e%d]  (log-scale)\n', lb(5), ub(5));
fprintf('-------------------------------------------\n');

%% ---- Run GCRA ----
tic;
[Best_score, Best_pos, CG_curve] = gcra(SearchAgents_no, Max_iteration, lb, ub, dim, fobj);
time_taken = toc;

%% ---- Results ----
Q1_opt = 10^Best_pos(5);

fprintf('\n=== OPTIMAL PARAMETERS FOUND (NONLINEAR) ===\n');
fprintf('Composite cost J = %g\n', Best_score);
fprintf('Time elapsed    = %.1f s\n\n', time_taken);

fprintf('--- Copy these values into get_nonlinear_params.m ---\n');
fprintf('p.Kp  = %f;\n', Best_pos(1));
fprintf('p.Ki  = %f;\n', Best_pos(2));
fprintf('p.Kd  = %f;\n', Best_pos(3));
fprintf('p.tau = %f;  %% IEID LPF time constant\n', Best_pos(4));
fprintf('Q     = diag([%g, 1e2, 10, 10, 1]);  %% Q1 = 10^%f\n', Q1_opt, Best_pos(5));
fprintf('p.L_ieid = lqe(p.A, eye(5), p.C, Q, 1);\n');
fprintf('-----------------------------------------------------\n');

%% ---- Convergence Plot ----
figure('Name','GCRA Convergence (Nonlinear)','Color','w');
semilogy(CG_curve, 'LineWidth', 2);
grid on;
xlabel('Iteration');
ylabel('Best Score J');
title('GCRA Convergence Curve for Nonlinear IEID');
