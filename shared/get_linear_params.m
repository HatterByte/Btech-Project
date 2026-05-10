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

%% IEID observer gain (from original params_pwr.m)
L = [1.957470876325828e+03;
     0.156014381434970;
     87.431287944182941;
     7.555735482558515;
     1.000000000000067e+02];

p.L_obs  = L;
p.A_obs  = A_n - L*C;
p.B_obs  = [B_n, L];    % inputs: [u_f; y]
p.Bnp    = pinv(B_n);
p.tau    = 0.02;

%% PID gains (outer loop)
p.Kp = 8.89;
p.Ki = 5.26;
p.Kd = 0.0521;
p.N  = 100;

end
