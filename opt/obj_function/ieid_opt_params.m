function p = ieid_opt_params(Kp, Ki, Kd, tau, Q1)
%IEID_OPT_PARAMS  Build IEID parameter struct for a given set of tunable values.
%
%   Returns empty [] if LQE computation fails (unstable observer placement).

% Base plant parameters
base = get_linear_params();

p = base;
p.Kp  = Kp;
p.Ki  = Ki;
p.Kd  = Kd;
p.N   = 100;     % Derivative filter — fixed, not tuned
p.tau = tau;

% Recompute LQE observer gain with tunable Q1
Q = diag([Q1, 1e2, 10, 10, 1]);
R = 1;

try
    L = lqe(p.A_real, eye(5), p.C, Q, R);
    p.L_obs = L;
    p.A_obs = p.A_real - L * p.C;
    p.B_obs = [p.B_real, L];
    p.Bnp   = pinv(p.B_real);
catch
    p = [];  % LQE failed — infeasible point
    return;
end

end
