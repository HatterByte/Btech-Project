function u = adrc_control_law(r, x_hat, p)
%ADRC_CONTROL_LAW Paper-matched ADRC/CESO feedback law.
%
%   x_hat = [estimated output; estimated output derivative; total disturbance]

u = (p.wc^2 * r ...
     + (-p.wc^2 + p.a0) * x_hat(1) ...
     + (-2*p.wc + p.a1) * x_hat(2) ...
     - x_hat(3)) / p.b0;
end
