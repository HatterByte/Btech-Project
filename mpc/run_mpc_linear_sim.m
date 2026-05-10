function [t_out, y_out, u_out, x_out] = run_mpc_linear_sim(p, mpcobj, tspan, x0)
%RUN_MPC_LINEAR_SIM  Simulates the linearized PWR plant with MPC control.
%   Uses a hybrid discrete-continuous approach: MPC updates at Ts, 
%   plant integrates continuously in between.

Ts = mpcobj.Ts;
t_start = tspan(1);
t_end   = tspan(2);
t_steps = t_start:Ts:t_end;

% Initialize
x = x0;
t_out = [];
y_out = [];
u_out = [];
x_out = [];

% MPC State
ms = mpcstate(mpcobj);

% Simulation Loop
fprintf('Running MPC simulation (discrete-continuous hybrid)...\n');
for k = 1:length(t_steps)-1
    t_now = t_steps(k);
    t_next = t_steps(k+1);
    
    % 1. Get current output (neutron density deviation)
    y = p.C * x;
    
    % 2. Get current reference
    r = 0.1 * (t_now >= 20);
    
    % 3. Call MPC to get control input u
    u = mpcmove(mpcobj, ms, y, r);
    
    % 4. Integrate plant for Ts with constant u
    % Note: disturbance is time-varying, so we still pass it to the ODE
    [t_ode, x_ode] = ode15s(@(t,x) plant_ode_with_dist(t, x, u, p), ...
                           [t_now t_next], x, odeset('RelTol', 1e-6));
    
    % 5. Store results
    t_out = [t_out; t_ode];
    x_out = [x_out; x_ode];
    u_out = [u_out; u * ones(length(t_ode), 1)];
    y_out = [y_out; (p.C * x_ode')'];
    
    % 6. Update plant state for next iteration
    x = x_ode(end, :)';
end

end

function dx = plant_ode_with_dist(t, x, u, p)
    % Plant: ẋ = A_real*x + B_real*u + B_d*d
    d = disturbance_fcn(t);
    dx = p.A_real * x + p.B_real * u + p.B_d * d;
end
