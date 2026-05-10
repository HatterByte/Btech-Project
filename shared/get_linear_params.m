function p = get_linear_params()
%GET_LINEAR_PARAMS Returns parameters for the OLD linearized PWR plant.
%
%   This is a struct version of the original params_pwr.m, for use with
%   the PID-only uncontrolled baseline simulation on the old linear plant.

%% Physical parameters
beta   = 0.006019;
Lambda = 2e-5;
lambda = 0.15;
ff     = 0.92;
P0     = 2500;
Gr     = 0.0145;
n0     = 0.5;

%% Derived
af    = (n0 - 4.24) / 1e5;
ac    = (-4*n0 - 17.3) / 1e5;
Mf    = 26.3;
Mc    = (160/9)*n0 + 54.002;
omega = (5/3)*n0 + 4.9333;
M     = 28*n0 + 74;

%% State-space matrices (5×5, deviation form around n0=0.5)
A_n = [ -beta/Lambda,      beta/Lambda,  (af*n0)/Lambda,  (ac*n0)/(2*Lambda),  n0/Lambda;
         lambda,           -lambda,       0,               0,                   0;
         (ff*P0)/Mf,        0,           -omega/Mf,        omega/(2*Mf),        0;
         ((1-ff)*P0)/Mf,    0,            omega/Mc,       -(2*M+omega)/(2*Mc),  0;
         0,                 0,            0,               0,                   0 ];

B_n = [0; 0; 0; 0; Gr];
C   = [1 0 0 0 0];
D   = 0;

B_d = B_n;     % disturbance enters through same channel as control

%% No model uncertainty
p.A_real = A_n;
p.B_real = B_n;
p.B_d    = B_d;
p.C      = C;
p.D      = D;

%% IEID observer gain (Optimal from GCRA)
Q = diag([1.11153, 1e2, 10, 10, 1]);  % Q1 = 10^0.0459
L = lqe(p.A_real, eye(5), p.C, Q, 1);

p.L_obs  = L;
p.A_obs  = p.A_real - L*p.C;
p.B_obs  = [p.B_real, L];    % inputs: [u_f; y]
p.Bnp    = pinv(p.B_real);
p.tau    = 0.004165;  % IEID LPF time constant

%% PID gains (outer loop)
p.Kp = 4.152582;
p.Ki = 2.708109;
p.Kd = 0.009036;
p.N  = 100;

end
