function [Score, Position, Convergence] = gcra(Search_Agents, Max_iterations, Lower_bound, Upper_bound, dimension, objective)
%GCRA  Golden Crane Rat Algorithm (modified for better exploration).
%
%   Improvements over original:
%     1. GR_r floor: minimum search radius prevents complete stagnation
%     2. Diversity restart: if no improvement for 'stall_limit' iters,
%        reinitialise the worst 30% of agents randomly near the best position
%     3. Alpha_score updated only after full iteration (not mid-loop),
%        which prevents bias toward early-iteration positions

Position    = zeros(1, dimension);
Score       = inf;
Gcanerats   = init(Search_Agents, dimension, Upper_bound, Lower_bound);
Convergence = zeros(1, Max_iterations);
l           = 1;
Alpha_pos   = Position;
Alpha_score = Score;

%% --- Initial population evaluation ---
for i = 1:size(Gcanerats, 1)
    Gcanerats(i,:) = min(max(Gcanerats(i,:), Lower_bound), Upper_bound);
    fitness = objective(Gcanerats(i,:));
    if fitness < Score
        Score       = fitness;
        Position    = Gcanerats(i,:);
        Alpha_pos   = Position;
        Alpha_score = Score;
    end
end

fprintf('  Init      | Best J = %.6f | Pos = [%s]\n', Score, num2str(Position, '%.4f  '));

%% --- Parameters ---
stall_limit = 5;       % trigger restart after this many non-improving iters
stall_count = 0;       % current stall counter
r_min_frac  = 0.02;    % GR_r floor = 2% of initial Alpha_score
r_min       = r_min_frac * Alpha_score;

%% --- Main loop ---
while l < Max_iterations + 1

    GR_m   = randperm(Search_Agents - 1, 1);
    GR_rho = 0.5;

    % Search radius: linear decay with a minimum floor to prevent stagnation
    GR_r = max(Alpha_score * (1 - l/Max_iterations), r_min);

    GR_mu    = floor(3 * rand + 1);     % uniform in {1,2,3}
    GR_c     = rand;
    GR_alpha = 2*GR_r*rand - GR_r;
    GR_beta  = 2*GR_r*GR_mu - GR_r;

    %% --- Attraction step: pull all agents toward Alpha ---
    for i = 1:size(Gcanerats, 1)
        for j = 1:size(Gcanerats, 2)
            Gcanerats(i,j) = (Gcanerats(i,j) + Alpha_pos(j)) / 2;
            Gcanerats(i,j) = min(max(Gcanerats(i,j), Lower_bound(j)), Upper_bound(j));
        end
    end

    %% --- Update step ---
    iter_best_score = Score;
    for i = 1:size(Gcanerats, 1)
        for j = 1:size(Gcanerats, 2)
            if rand < GR_rho
                % Exploitation branch A: move toward Alpha
                Gcanerats(i,j) = Gcanerats(i,j) + GR_c*(Alpha_pos(j) - GR_r*Gcanerats(i,j));
                Gcanerats(i,j) = min(max(Gcanerats(i,j), Lower_bound(j)), Upper_bound(j));

                fitness = objective(Gcanerats(i,:));
                if fitness < iter_best_score
                    iter_best_score = fitness;
                    Position = Gcanerats(i,:);
                else
                    % Fallback A: exploration perturbation
                    Gcanerats(i,j) = Gcanerats(i,j) + GR_c*(Gcanerats(i,j) - GR_alpha*Alpha_pos(j));
                    Gcanerats(i,j) = min(max(Gcanerats(i,j), Lower_bound(j)), Upper_bound(j));

                    fitness = objective(Gcanerats(i,:));
                    if fitness < iter_best_score
                        iter_best_score = fitness;
                        Position = Gcanerats(i,:);
                    end
                end
            else
                % Exploitation branch B: interact with random agent
                Gcanerats(i,j) = Gcanerats(i,j) + GR_c*(Alpha_pos(j) - GR_mu*Gcanerats(GR_m,j));
                Gcanerats(i,j) = min(max(Gcanerats(i,j), Lower_bound(j)), Upper_bound(j));

                fitness = objective(Gcanerats(i,:));
                if fitness < iter_best_score
                    iter_best_score = fitness;
                    Position = Gcanerats(i,:);
                else
                    % Fallback B
                    Gcanerats(i,j) = Gcanerats(i,j) + GR_c*(Gcanerats(GR_m,j) - GR_beta*Alpha_pos(j));
                    Gcanerats(i,j) = min(max(Gcanerats(i,j), Lower_bound(j)), Upper_bound(j));

                    fitness = objective(Gcanerats(i,:));
                    if fitness < iter_best_score
                        iter_best_score = fitness;
                        Position = Gcanerats(i,:);
                    end
                end
            end
        end
    end

    %% --- Update global best ---
    if iter_best_score < Score
        Score       = iter_best_score;
        Alpha_pos   = Position;
        Alpha_score = Score;
        stall_count = 0;
    else
        stall_count = stall_count + 1;
    end

    %% --- Diversity restart if stalled ---
    if stall_count >= stall_limit
        n_restart = max(1, floor(0.3 * Search_Agents));  % restart worst 30%
        % Reinitialise randomly, biased near best position (±20% of range)
        for i = 1:n_restart
            for j = 1:dimension
                span  = 0.2 * (Upper_bound(j) - Lower_bound(j));
                new_val = Alpha_pos(j) + span*(2*rand-1);
                Gcanerats(i,j) = min(max(new_val, Lower_bound(j)), Upper_bound(j));
            end
        end
        stall_count = 0;
    end

    %% --- Log iteration ---
    l           = l + 1;
    Convergence(l) = Score;
    fprintf('  Iter %3d/%d | Best J = %.6f | GR_r = %.4f | Stall = %d | Pos = [%s]\n', ...
            l-1, Max_iterations, Score, GR_r, stall_count, num2str(Position, '%.4f  '));
end
end

%% --- Initialisation helper ---
function Pos = init(SearchAgents, dimension, upperbound, lowerbound)
Boundary = size(upperbound, 2);
if Boundary == 1
    Pos = rand(SearchAgents, dimension) .* (upperbound - lowerbound) + lowerbound;
end
if Boundary > 1
    for i = 1:dimension
        Pos(:,i) = rand(SearchAgents,1) .* (upperbound(i) - lowerbound(i)) + lowerbound(i);
    end
end
end