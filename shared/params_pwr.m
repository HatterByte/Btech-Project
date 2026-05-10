%% ============================================================
%%  PWR SISO MODEL 
%% ============================================================

%% ---------- Parameters ----------
beta = 0.006019;
Lambda = 2e-5;
lambda = 0.15;
ff = 0.92;
P0 = 2500;   

Gr = 0.0145;    

n0 = 0.5;
% n0 = 1;
af = (n0-4.24)/100000;
ac = (-4*n0 - 17.3)/100000;
Mf = 26.3;
Mc = (160/9*n0 + 54.002);
omega = (5/3*n0 + 4.9333);
M = (28*n0 + 74);


%% ---------- State-space matrices ----------
% A_n is 5×5
A_n = [(-beta/Lambda) (beta/Lambda) (af*n0/Lambda) (ac*n0/(2*Lambda)) n0/Lambda;
        lambda -lambda 0 0 0;
        (ff*P0/Mf) 0 (-omega/Mf) (omega/(2*Mf)) 0;
        ((1-ff)*P0/Mf) 0 (omega/Mc) (-(2*M + omega)/(2*Mc)) 0;
        0 0 0 0 0];

% B is 5×1
B_n = [ 0; 0; 0; 0; Gr ];

% Output is neutron density 
C = [1 0 0 0 0];

D = 0;
  

%% ---------- Disturbance channel ----------
B_d = B_n;   
%% ---------- uncertainty ----------
DeltaA = zeros(5);
DeltaB = zeros(5,1);

%% ---------- Real plant (for simulation) ----------
A_real = A_n + DeltaA;
B_real = B_n + DeltaB;

eigA = eig(A_n);
ContRank = rank(ctrb(A_n,B_n));
ObsRank  = rank(obsv(A_n,C));
% disp(eigA);

%% ============================================================
%%  Observer gain 
%% ============================================================






% Q_scales = logspace(0,8,20);
% R_scales = logspace(-4,2,20);
% 
% best = inf;
% bestL = [];

% for q = Q_scales
%     for r = R_scales
%         Q = diag([q q/100 10 10 1]);
%         R = r;
% 
%         L_try = lqe(A_n, eye(5), C, Q, R);
%         poles = eig(A_n - L_try*C);
% 
%         % Reject unstable observers
%         if any(real(poles) > 0)
%             continue;
%         end
% 
%         % Evaluate performance (tracking + disturbance)
%         J = evalObserverPerformance(L_try); % You write this small function
% 
%         if J < best
%             best = J;
%             bestL = L_try;
%         end
%     end
% end


L = [1.957470876325828e+03,
0.156014381434970,
87.431287944182941,
7.555735482558515,
1.000000000000067e+02];






% LQE alternative
% Qw = diag([1e6, 1e3, 10, 10, 0.1]);   % process noise (tune)
% Rv = 1;                               % measurement noise (tune) 
% L_kf = lqe(A_n, eye(size(A_n)), C, Qw, Rv);  % returns estimator gain
% % fprintf('norm(L_kf) = %.3e\n', norm(L_kf));
% eig_obs = eig(A_n - L_kf*C);
% L = L_kf;
%% ============================================================
%%  IEID estimator
%% ============================================================
Bnp = pinv(B_n);

tau = 0.02;
numF = [1];
denF = [tau^2  1.41*tau  1];

%% ============================================================
%%  PID gains (outer loop)
%% ============================================================
% Kp = 10;
% Ki = 4.7243;
% Kd = 0.28253;


Kp = 8.89;
Ki = 5.26;
Kd = 0.0521;
%% ============================================================
%%  Combined matrices
%% ============================================================
B_total = [B_real  B_d];
D_total = zeros(size(C,1),2);

%% ---------- Observer system ----------
A_obs = A_n - L*C;
B_obs = [B_n  L];   % Inputs = [u ; y]
C_obs = eye(5);
D_obs = zeros(5,2);



