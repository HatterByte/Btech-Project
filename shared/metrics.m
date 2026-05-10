function m = metrics(t, y, r, u)
%METRICS Computes performance metrics for the control system.
%
%   INPUTS:
%     t : Time vector
%     y : Output signal
%     r : Reference signal
%     u : Control effort

e = r - y;
m.rms_error = sqrt(mean(e.^2));
m.control_energy = sum(u.^2) * mean(diff(t));

% Rise time (10% to 90% of the step change)
% Step occurs at t=20, final target is 0.1
target = 0.1;
t_step = 20;
idx_step = find(t >= t_step, 1);

t_after = t(idx_step:end);
y_after = y(idx_step:end);

y_start = 0;
y_final = target;
delta_y = y_final - y_start;

idx10 = find(y_after >= y_start + 0.1*delta_y, 1);
idx90 = find(y_after >= y_start + 0.9*delta_y, 1);

if ~isempty(idx10) && ~isempty(idx90)
    m.rise_time = t_after(idx90) - t_after(idx10);
else
    m.rise_time = NaN;
end

% Settling time (2% band)
band = 0.02 * delta_y;
idx_settle = find(abs(y_after - y_final) > band, 1, 'last');
if ~isempty(idx_settle) && idx_settle < length(t_after)
    m.settling_time = t_after(idx_settle) - t_step;
else
    m.settling_time = 0; % Already settled or never left
end

% Overshoot
m.overshoot = (max(y_after) - y_final) / delta_y * 100;
if m.overshoot < 0, m.overshoot = 0; end

% Steady-state error (mean of last 5 seconds)
idx_ss = t > (t(end) - 5);
m.sse = mean(abs(r(idx_ss) - y(idx_ss)));

end
