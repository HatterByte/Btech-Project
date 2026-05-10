function mpcobj = setup_mpc(p)
%SETUP_MPC  Initializes the MPC controller for the linearized PWR plant.
%
%   INPUTS:
%     p - parameter struct with A_real, B_real, C, D
%
%   OUTPUTS:
%     mpcobj - MATLAB MPC object

% 1. Define the continuous plant model
sys_c = ss(p.A_real, p.B_real, p.C, p.D);

% 2. Discretize for MPC
Ts = 0.1;   % 100ms sampling time
sys_d = c2d(sys_c, Ts);

% 3. Create MPC object
% Prediction horizon Hp = 20 steps (2.0 seconds)
% Control horizon Hc = 3 steps (0.3 seconds)
Hp = 20;
Hc = 3;

mpcobj = mpc(sys_d, Ts, Hp, Hc);

% 4. Weights
% We want tight power tracking (y1 = dn)
mpcobj.Weights.OutputVariables = 1.0;
% We want to penalize excessive rod movement rate (u)
mpcobj.Weights.ManipulatedVariables = 0; % No penalty on absolute position deviation
mpcobj.Weights.ManipulatedVariablesRate = 0.05;

% 5. Constraints (from rod velocity limits in ADRC paper)
mpcobj.MV.Min = -0.01;
mpcobj.MV.Max =  0.01;

% 6. Initial state
% The plant is in deviation form, so equilibrium is 0.
% Note: We don't need to set this here, mpcmove handles it.

end
