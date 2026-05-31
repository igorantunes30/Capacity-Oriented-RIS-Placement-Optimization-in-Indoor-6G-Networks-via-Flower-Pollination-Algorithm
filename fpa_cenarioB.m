function [best, fbest, convergencia, sc_best] = fpa_cenarioB(para, N_elementos, pos_Rx, cfg, csv_path)
% FPA para o CENARIO B. A RIS percorre as paredes laterais Sul (y=0) e
% Norte (y=50) -- ambas Scenario 1. Avaliacao via calculamediavazao_cenarioB.

    if nargin < 4 || isempty(cfg); cfg = cenarioB_config(); end

    n = para(1);
    p = para(2);
    d = 3;
    N_iter = 300;
    gamma_step = 0.1;
    beta = 1.5;

    Lb = cfg.Lb; Ub = cfg.Ub;

    Sol = zeros(n, d);
    Scen = zeros(n, 1);
    for i = 1:n
        cand0 = Lb + rand(1,d) .* (Ub - Lb);
        [Sol(i,:), Scen(i)] = snap_lateral(cand0, Lb, Ub);
    end

    Fitness = zeros(n, 1);
    parfor i = 1:n
        Fitness(i) = calculamediavazao_cenarioB(Sol(i,:), Scen(i), N_elementos, pos_Rx, cfg);
    end
    [fbest, I] = max(Fitness);
    best = Sol(I,:);
    sc_best = Scen(I);

    convergencia = zeros(N_iter, 1);
    sigma = (gamma(1+beta)*sin(pi*beta/2) / (gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);

    for t = 1:N_iter
        S = zeros(n, d);
        Ssc = zeros(n, 1);
        for i = 1:n
            if rand > p
                u = randn(1, d) * sigma;
                v = randn(1, d);
                Lf = u ./ abs(v).^(1/beta);
                cand = Sol(i,:) + gamma_step * Lf .* (Sol(i,:) - best);
            else
                epsilon = rand;
                JK = randperm(n);
                cand = Sol(i,:) + epsilon * (Sol(JK(1),:) - Sol(JK(2),:));
            end
            [S(i,:), Ssc(i)] = snap_lateral(cand, Lb, Ub);
        end

        Fnew = zeros(n, 1);
        parfor i = 1:n
            Fnew(i) = calculamediavazao_cenarioB(S(i,:), Ssc(i), N_elementos, pos_Rx, cfg);
        end

        for i = 1:n
            if Fnew(i) >= Fitness(i)
                Sol(i,:) = S(i,:);
                Fitness(i) = Fnew(i);
            end
            if Fnew(i) > fbest
                best = S(i,:);
                fbest = Fnew(i);
                sc_best = Ssc(i);
            end
        end
        convergencia(t) = fbest;

        if nargin >= 5 && ~isempty(csv_path)
            fid = fopen(csv_path, 'a');
            if fid > 0
                fprintf(fid, '%d,%d,%.6f,%.3f,%.3f,%.3f\n', ...
                        N_elementos, t, fbest, best(1), best(2), best(3));
                fclose(fid);
            end
        end

        if t == 1 || mod(t, 25) == 0
            fprintf('   iter %4d/%d | melhor fitness = %.4f | pos [%.1f, %.1f, %.1f]\n', ...
                    t, N_iter, fbest, best(1), best(2), best(3));
        end
    end
end

% -------------------------------------------------------------------------
function [s, sc] = snap_lateral(s, Lb, Ub)
    % Projeta na parede lateral mais proxima (Sul y=0 ou Norte y=50).
    s = max(s, Lb);
    s = min(s, Ub);
    if (s(2) - 0) <= (Ub(2) - s(2))
        s(2) = 0;  sc = 1;   % Sul
    else
        s(2) = Ub(2); sc = 1; % Norte
    end
end
