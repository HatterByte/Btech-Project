function p = get_paper_params()
%GET_PAPER_PARAMS Parameters exactly as defined in the ADRC paper.
%
%   Reference: "Active disturbance rejection control of pressurized water reactor"
%   S. Ahmad et al., Annals of Nuclear Energy, 189, 109845, 2023.
%
%   Parameters confirmed from paper's Simulink source code (images):
%     tot_react.m  : Tf0=650, Tl0=314, Te0=290, af=(n-4.24)*1e-5, ac=(-4n-17.3)*1e-5
%     npp_nl.m     : B1=0.8*B, lmb1=0.8*lmb, l=0.0001
%     fb_lin.m     : Feedback linearization, Gr1=1.2*Gr

%% --- Physical Constants (from npp_nl) ---
p.beta   = 0.006019;       % B  — delayed neutron fraction
p.Lambda = 1e-4;           % l  — prompt neutron lifetime (s)
p.lambda = 0.15;           % lmb — precursor decay constant (s^-1)
p.ff     = 0.92;           % ff — fraction of fission heat to fuel
p.P0     = 2500;           % P  — nominal thermal power (MW)
p.muf    = 26.3;           % uf — fuel thermal capacity (MWs/°C)
p.Gr     = 0.0145;         % Gr — rod worth rate (δk/k/s per unit u)
p.Gr1    = 1.2 * 0.0145;   % Gr1 = 1.2*Gr (used in fb_lin)

% Uncertain model parameters (as in paper)
p.B1   = 0.8 * p.beta;     % B1  = 0.8*B
p.lmb1 = 0.8 * p.lambda;   % lmb1 = 0.8*lmb

%% --- Nominal Operating Point: 100% power ---
p.n0 = 1.0;
p.c0 = 1.0;

%% --- Reference temperatures (from tot_react.m — confirmed from paper image) ---
p.Tf0 = 650;    % Tf0 — fuel temperature reference (°C)
p.Tl0 = 314;    % Tl0 — coolant temperature reference (°C)
p.Te0 = 290;    % Te0 — inlet temperature reference (°C)

% These are used in the reactivity formula and also as initial conditions.
% The plant equilibrium at n=1.0 solves to Tf≈650.7°C, Tl≈314.5°C —
% matching these values to within rounding (< 0.5°C error).
p.Tf_ref = p.Tf0;
p.Tl_ref = p.Tl0;
p.Te_ref = p.Te0;

%% --- ADRC / ESO Parameters (confirmed from paper's Simulink ESO block) ---
% From ESO code:
%   bo=150; k=1; wo=55*k; wc=1;
%   L = [3*wo; 3*wo^2; wo^3]  (standard bandwidth parameterization)
%   u = sat( (wc^2*(ref-z1) - 2*wc*z2 - z3) / bo )
p.b0 = 150;          % bo — confirmed from paper ESO code (bo=150)
p.wo = 65;           % observer bandwidth — increased from 55 (Suggestion 1: tighter tracking)
p.wc = 1;            % controller bandwidth — confirmed from paper (wc=1)
p.k1 = p.wc^2;       % = 1
p.k2 = 2 * p.wc;     % = 2

%% --- Rod velocity limits ---
p.u_max =  0.01;
p.u_min = -0.01;

end

