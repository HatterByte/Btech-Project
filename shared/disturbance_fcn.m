function d = disturbance_fcn(t)
%DISTURBANCE_FCN  External disturbance applied to the PWR plant input.
%
%   d = 0.002 * ( sin(5t) + 0.5*tanh(4*(t-20)) + cos(3t) )
%
%   Components:
%     sin(5t)          - sinusoidal oscillation
%     0.5*tanh(4(t-20))- smooth step-like transition at t = 20 s
%     cos(3t)          - additional oscillatory disturbance

    d = 0.002 * ( sin(5*t) + 0.5*tanh(4*(t - 20)) + cos(3*t) );
end
