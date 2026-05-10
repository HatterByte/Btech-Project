function dz = run_closed_loop(t, z, p)
%RUN_CLOSED_LOOP  ODE right-hand-side for the full IEID-PWR closed-loop.
%
%   Augmented state vector z (14 x 1):
%     z(1:5)   - plant states  x  = [dn, dc, dTf, dTl, drho_r]'
%     z(6:10)  - observer states  x_hat
%     z(11)    - PID integrator state  (integral of error)
%     z(12)    - PID derivative filter state  (1st-order, N = 100)
%     z(13)    - IEID LPF state 1  (filter output  = d_tilde)
%     z(14)    - IEID LPF state 2  (filter derivative state)
%
%   Signal flow (per time step):
%     r_dev  = r(t) - 0.5
%     y      = C * x_plant                         (plant output)
%     e      = r_dev - y                           (tracking error)
%     u_f    = Kp*e + Ki*x_int + Kd*N*(e - x_d)   (PID, filtered deriv)
%     d_tilde= z(13)                               (LPF output from prev step)
%     u      = u_f - d_tilde                       (compensated input)
%     d      = disturbance_fcn(t)                  (external disturbance)
%     d_hat  = IEID_estimator(u_f, u, x_hat, y)   (raw estimate)
%     -> LPF: d_hat filtered -> d_tilde (states z(13:14))
%
%   PARAMETER STRUCT p (set by main_ieid_pwr.m):
%     p.A_real, p.B_real, p.B_d
%     p.A_obs,  p.B_obs
%     p.C,      p.Bnp,  p.L
%     p.Kp, p.Ki, p.Kd, p.N
%     p.tau

    %% ---- Unpack states ------------------------------------------------
    x_plant = z(1:5);       % plant states
    x_hat   = z(6:10);      % observer states
    x_int   = z(11);        % PID integral
    x_d     = z(12);        % PID derivative filter
    xf1     = z(13);        % LPF output  (= d_tilde)
    xf2     = z(14);        % LPF rate state

    %% ---- Reference signal ---------------------------------------------
    % Step: r = 0.5 for t < 20, r = 0.6 for t >= 20
    r     = 0.5 + 0.1 * (t >= 20);
    r_dev = r - 0.5;        % deviation from nominal operating point

    %% ---- Plant output -------------------------------------------------
    y = p.C * x_plant;      % scalar: neutron density deviation

    %% ---- Tracking error -----------------------------------------------
    e = r_dev - y;

    %% ---- PID controller (Simulink parallel form, N = 100 filter) ------
    % Derivative filter:  H_d(s) = N*s / (s + N)
    % State x_d: ẋ_d = -N*x_d + N*e
    % Derivative term   = Kd * N * (e - x_d)
    deriv_term = p.Kd * p.N * (e - x_d);
    u_f = p.Kp * e + p.Ki * x_int + deriv_term;

    %% ---- IEID compensation -------------------------------------------
    d_tilde = xf1;           % current filtered disturbance estimate
    u       = u_f - d_tilde; % compensated plant input

    %% ---- External disturbance ----------------------------------------
    d = disturbance_fcn(t);

    %% ---- IEID raw estimator ------------------------------------------
    d_hat = ieid_estimator_fcn(u_f, u, x_hat, y, p.L, p.C, p.Bnp);

    %% ---- Compute derivatives -----------------------------------------

    % 1) Plant:  ẋ = A_real*x + B_real*u + B_d*d
    dx_plant = p.A_real * x_plant + p.B_real * u + p.B_d * d;

    % 2) Observer:  ẋ_hat = A_obs*x_hat + B_obs*[u_f; y]
    %               where A_obs = A_n - L*C,  B_obs = [B_n, L]
    dx_hat = p.A_obs * x_hat + p.B_obs * [u_f; y];

    % 3) PID integrator:  ẋ_int = e
    dx_int = e;

    % 4) PID derivative filter:  ẋ_d = -N*x_d + N*e
    dx_d = -p.N * x_d + p.N * e;

    % 5) IEID 2nd-order LPF:  τ²ÿ + 1.41τẏ + y = d_hat
    %    States: xf1 = y_filt (= d_tilde),  xf2 = ẏ_filt
    dxf1 = xf2;
    dxf2 = (d_hat - xf1 - 1.41 * p.tau * xf2) / (p.tau^2);

    %% ---- Assemble derivative vector ----------------------------------
    dz = [dx_plant; dx_hat; dx_int; dx_d; dxf1; dxf2];
end
