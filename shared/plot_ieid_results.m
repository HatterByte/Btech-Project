function plot_ieid_results(t, sig)
%PLOT_IEID_RESULTS  Generate publication-quality plots for IEID-PWR simulation.
%
%   INPUTS
%     t   : Time vector (Nx1)
%     sig : Struct with logged signals:
%           sig.r        - reference (Nx1)
%           sig.r_dev    - reference deviation (Nx1)
%           sig.y        - plant output (Nx1)
%           sig.e        - tracking error (Nx1)
%           sig.u_f      - PID output (Nx1)
%           sig.u        - compensated plant input (Nx1)
%           sig.d        - actual external disturbance (Nx1)
%           sig.d_hat    - raw IEID estimate (Nx1)
%           sig.d_tilde  - filtered IEID estimate (Nx1)

    grey  = [0.45 0.45 0.45];
    blue  = [0.08 0.40 0.75];
    red   = [0.85 0.15 0.15];
    green = [0.10 0.65 0.30];
    orng  = [0.90 0.45 0.00];

    lw = 1.8;
    fs = 12;

    %% ============================================================
    %% Figure 1 — Power tracking (r_dev vs y)
    %% ============================================================
    figure('Name','Power Tracking','NumberTitle','off', ...
           'Position',[100 500 820 340]);

    plot(t, sig.r_dev, '--', 'Color', grey,  'LineWidth', lw, ...
         'DisplayName','Reference r_{dev}'); hold on;
    plot(t, sig.y,     '-',  'Color', blue,  'LineWidth', lw, ...
         'DisplayName','Plant output y');
    hold off;

    xlabel('Time (s)', 'FontSize', fs);
    ylabel('Neutron density deviation', 'FontSize', fs);
    title('Reactor Power Tracking — IEID Control', 'FontSize', fs+1);
    legend('Location','southeast','FontSize',fs-1);
    grid on; xlim([t(1) t(end)]);
    set(gca,'FontSize',fs);

    %% ============================================================
    %% Figure 2 — Disturbance estimation
    %% ============================================================
    figure('Name','Disturbance Estimation','NumberTitle','off', ...
           'Position',[100 100 820 340]);

    plot(t, sig.d,       '-',  'Color', red,   'LineWidth', lw, ...
         'DisplayName','Actual disturbance d(t)'); hold on;
    plot(t, sig.d_hat,   '--', 'Color', orng,  'LineWidth', lw-0.3, ...
         'DisplayName','Raw estimate \hat{d}(t)');
    plot(t, sig.d_tilde, '-',  'Color', green, 'LineWidth', lw, ...
         'DisplayName','Filtered estimate \tilde{d}(t)');
    hold off;

    xlabel('Time (s)', 'FontSize', fs);
    ylabel('Disturbance (equivalent input)', 'FontSize', fs);
    title('IEID Disturbance Estimation', 'FontSize', fs+1);
    legend('Location','best','FontSize',fs-1);
    grid on; xlim([t(1) t(end)]);
    set(gca,'FontSize',fs);

    %% ============================================================
    %% Figure 3 — Control inputs
    %% ============================================================
    figure('Name','Control Inputs','NumberTitle','off', ...
           'Position',[940 500 820 340]);

    plot(t, sig.u_f, '-',  'Color', blue,  'LineWidth', lw, ...
         'DisplayName','PID output u_f'); hold on;
    plot(t, sig.u,   '--', 'Color', red,   'LineWidth', lw, ...
         'DisplayName','Compensated input u');
    hold off;

    xlabel('Time (s)', 'FontSize', fs);
    ylabel('Control signal (rod reactivity rate)', 'FontSize', fs);
    title('Control Inputs — PID vs Compensated', 'FontSize', fs+1);
    legend('Location','best','FontSize',fs-1);
    grid on; xlim([t(1) t(end)]);
    set(gca,'FontSize',fs);

    %% ============================================================
    %% Figure 4 — Tracking error
    %% ============================================================
    figure('Name','Tracking Error','NumberTitle','off', ...
           'Position',[940 100 820 340]);

    plot(t, sig.e, '-', 'Color', red, 'LineWidth', lw);
    xlabel('Time (s)', 'FontSize', fs);
    ylabel('Error e(t) = r_{dev} - y', 'FontSize', fs);
    title('Tracking Error', 'FontSize', fs+1);
    grid on; xlim([t(1) t(end)]);
    set(gca,'FontSize',fs);
    yline(0,'--','Color',grey,'LineWidth',1);

end
