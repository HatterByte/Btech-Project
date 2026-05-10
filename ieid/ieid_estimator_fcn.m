function d_hat = ieid_estimator_fcn(u_f, u, x_hat, y, L, C, Bnp)
%IEID_ESTIMATOR_FCN  Equivalent Input Disturbance (IEID) estimator.
%
%   Computes the raw disturbance estimate d_hat from:
%     1. Observer output residual  : y - C*x_hat
%     2. Existing compensation term: u_f - u  (= d_tilde fed back)
%
%   INPUTS
%     u_f   : PID controller output (scalar)
%     u     : Actual plant input = u_f - d_tilde (scalar)
%     x_hat : Observer state estimate (5x1)
%     y     : Measured plant output (scalar)
%     L     : Observer gain (5x1)
%     C     : Output matrix (1x5)
%     Bnp   : Pseudo-inverse of B_n (1x5)
%
%   OUTPUT
%     d_hat : Raw disturbance estimate (scalar)

    % Step 1 — Observer predicted output
    y_hat = C * x_hat;

    % Step 2 — Output residual (real plant vs observer)
    residual = y - y_hat;

    % Step 3 — Convert residual to equivalent input disturbance space
    %          Bnp * L scales from output space back to input space
    scaled_error = Bnp * (L * residual);

    % Step 4 — Existing disturbance compensation already applied
    d_tilde_current = u_f - u;

    % Step 5 — Total raw disturbance estimate
    d_hat = scaled_error + d_tilde_current;
end
