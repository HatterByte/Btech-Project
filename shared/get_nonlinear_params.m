function p = get_nonlinear_params()
%GET_NONLINEAR_PARAMS Returns the parameters for the nonlinear PWR plant.
%   As specified in the ADRC paper/prompt.

p.beta   = 0.006019;
p.lambda = 0.15;
p.L      = 0.0001;
p.ff     = 0.92;
p.P      = 2500;
p.Gr     = 0.0145;

% Uncertain/modified kinetics parameters
p.B1   = 0.8 * p.beta;
p.lmb1 = 0.8 * p.lambda;

% Equilibrium point
p.n0 = 0.5;

% Derived parameters at equilibrium
p.M   = 28 * p.n0 + 74;
p.omg = (5/3) * p.n0 + 4.933;
p.uf  = 26.3;
p.uc  = (160/9) * p.n0 + 54.022;

% Feedback coefficients at equilibrium
p.af = (p.n0 - 4.24) * 1e-5;
p.ac = (-4 * p.n0 - 17.3) * 1e-5;

% Linearization at n0 = 0.5
% States: [dn; dc; dTf; dTl; drho_rod]
p.A = [ -p.B1/p.L, p.B1/p.L, (p.n0 * p.af)/p.L, (p.n0 * 0.5 * p.ac)/p.L, p.n0/p.L;
         p.lmb1,  -p.lmb1,   0,                  0,                      0;
         (p.ff * p.P)/p.uf, 0, -p.omg/p.uf, p.omg/(2*p.uf), 0;
         ((1-p.ff)*p.P)/p.uc, 0, p.omg/p.uc, -(2*p.M + p.omg)/(2*p.uc), 0;
         0, 0, 0, 0, 0 ];

p.B = [0; 0; 0; 0; p.Gr];
p.C = [1, 0, 0, 0, 0];
p.D = 0;

% Observer gain (Placeholder - will use one that works with new matrices)
% We'll use LQE or pole placement to find a stable L for this A matrix.
Q = diag([1e6, 1e2, 10, 10, 1]);
R = 1;
p.L_ieid = lqe(p.A, eye(5), p.C, Q, R);

% ADRC/CESO parameters from the paper Simulink setup.
p.b0  = 150;
p.k   = 1;
p.wo  = 55 * p.k;
p.wlg = 45 * p.k;
p.w1  = 20 * p.k;
p.w2  = 40 * p.k;
p.wc  = 1;

% Paper-specific reduced ADRC model.
p.a0 = -p.beta * p.lambda / p.L;
p.a1 =  p.beta / p.L;
p.Af = [0, 1, 0;
       -p.a0, -p.a1, 1;
        0, 0, 0];
p.Bf = [0; p.b0; 0];
p.Cf = [1, 0, 0];

% Standard ESO gain in the paper coordinates.
p.Lf = [3*p.wo - p.a1;
        3*p.wo^2 - p.a0 - p.a1*(3*p.wo - p.a1);
        p.wo^3];

% Two-level CESO gains.
p.L1 = [3*p.w1 - p.a1;
        3*p.w1^2 - p.a0 - p.a1*(3*p.w1 - p.a1);
        p.w1^3];
p.L2 = [3*p.w2 - p.a1;
        3*p.w2^2 - p.a0 - p.a1*(3*p.w2 - p.a1);
        p.w2^3];
p.Gamma = [0, 0, 0;
           0, 0, 1;
           0, 0, 0];

% Kept for compatibility with older scripts.
p.k1 = p.wc^2;
p.k2 = 2 * p.wc;

% PID Gains (for IEID outer loop)
p.Kp = 8.89;
p.Ki = 5.26;
p.Kd = 0.0521;
p.tau = 0.02;
p.N  = 100;  % Derivative filter coefficient (Simulink default)

% Control rod velocity limits
p.u_max =  0.01;
p.u_min = -0.01;

end
