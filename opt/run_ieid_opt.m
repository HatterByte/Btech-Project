%% RUN_IEID_OPT  —  Optimize IEID parameters using GCRA metaheuristic
%
%   Tunable parameters (5D search space):
%     x(1) = Kp          [1.0 , 25.0]
%     x(2) = Ki          [0.1 , 20.0]
%     x(3) = Kd          [0.0 ,  0.2]  — allows 0, optimizer decides
%     x(4) = tau         [0.001, 0.1]  — LPF time constant (prev run hit lb=0.005)
%     x(5) = log10(Q1)   [-1  ,  6 ]  — prev run hit lb=2 → extended down
%
%   Objective (ieid_objective.m) minimises composite cost:
%     J = ITAE(t≥20s) + overshoot penalty + SSE + rod saturation
%
%   GCRA improvements (in gcra.m):
%     - GR_r floor: never drops to zero, maintains exploration
%     - Diversity restart: reinitialises 30% of agents after stall_limit iters
%
%   Workflow:
%     1. Run this script once  →  GCRA finds optimal params
%     2. Copy printed values   →  paste into get_linear_params.m
%     3. Run main_linear.m     →  verify improvement

clear; clc;

%% ---- Paths ----
addpath('opt');
addpath('opt/core');
addpath('opt/obj_function');
addpath('shared');

%% ---- GCRA Settings ----
% Key insight from previous run:
%   - Most improvement happened in first 10-15 iterations
%   - Algorithm stalled from iter 8 onward (search radius → 0)
%   - Fix: more agents (better initial spread) + improved GCRA exploration
Search_Agents  = 50;
Max_iterations = 40;

%% ---- Search bounds ----
%        Kp     Ki      Kd      tau     log10(Q1)
lb = [  1.0,   0.1,   0.001,  0.003,      -1  ];
ub = [ 25.0,  20.0,   0.2,  0.100,       6  ];
dim = 5;

%   Bound rationale (from previous run analysis):
%     Kp/Ki : widened; previous optimal was Kp≈4.3, Ki≈2.7 (safe mid-range)
%     Kd    : [0, 0.2] — include 0 so optimizer can choose; prev hit lb=0.001
%     tau   : [0.001, 0.1] — prev run wanted tau<0.005; go lower
%     Q1    : [-1, 6] — prev run hit lower bound log10(Q1)=2; extend down

fitness = @ieid_objective;

%% ---- Run ----
fprintf('=== IEID Parameter Optimization (GCRA) ===\n');
fprintf('Agents: %d | Iterations: %d | Dimensions: %d\n', Search_Agents, Max_iterations, dim);
fprintf('Bounds:\n');
fprintf('  Kp  : [%.1f , %.1f]\n',  lb(1), ub(1));
fprintf('  Ki  : [%.1f , %.1f]\n',  lb(2), ub(2));
fprintf('  Kd  : [%.3f , %.3f]\n',  lb(3), ub(3));
fprintf('  tau : [%.4f, %.3f] s\n', lb(4), ub(4));
fprintf('  Q1  : [1e%.0f, 1e%.0f]  (log-scale)\n', lb(5), ub(5));
fprintf('-------------------------------------------\n');

tic;
[Score, Position, Convergence] = gcra(Search_Agents, Max_iterations, lb, ub, dim, fitness);
t_elapsed = toc;

%% ---- Results ----
Kp_opt  = Position(1);
Ki_opt  = Position(2);
Kd_opt  = Position(3);
tau_opt = Position(4);
Q1_opt  = 10^Position(5);

fprintf('\n=== OPTIMAL PARAMETERS FOUND ===\n');
fprintf('Composite cost J = %.6f\n', Score);
fprintf('Time elapsed    = %.1f s\n', t_elapsed);
fprintf('\n--- Copy these values into get_linear_params.m ---\n');
fprintf('p.Kp  = %.6f;\n', Kp_opt);
fprintf('p.Ki  = %.6f;\n', Ki_opt);
fprintf('p.Kd  = %.6f;\n', Kd_opt);
fprintf('p.tau = %.6f;  %% IEID LPF time constant\n', tau_opt);
fprintf('Q     = diag([%.6g, 1e2, 10, 10, 1]);  %% Q1 = 10^%.4f\n', Q1_opt, Position(5));
fprintf('p.L_ieid = lqe(p.A_real, eye(5), p.C, Q, 1);\n');
fprintf('--------------------------------------------------\n');

%% ---- Convergence plot ----
figure('Name','IEID Optimization — GCRA Convergence','Color','w','Position',[100 100 800 400]);
plot(Convergence(2:end), 'b-o', 'LineWidth', 1.8, 'MarkerSize', 5, 'DisplayName','Best J');
grid on;
xlabel('Iteration');
ylabel('Composite Cost J');
title('GCRA Convergence — IEID Parameter Optimization (5D)');
legend('Location','best');
