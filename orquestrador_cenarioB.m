% =========================================================================
% orquestrador_cenarioB.m
% Otimizacao do CENARIO B: sala assimetrica, paredes grossas e obstaculos
% bloqueando a visada direta dos 2 Tx. Usuarios aleatorios (semente fixa).
% RIS unica nas paredes laterais (Sul/Norte). FPA com 300 iteracoes.
% Gera curva de convergencia por N e salva resultado_cenarioB.mat.
% =========================================================================
clear; clc; close all;
clear functions;
rng(1);

% --- bootstrap de caminhos (cenario em pasta propria) -------------------
thisdir = fileparts(mfilename('fullpath'));
addpath(thisdir);             % funcoes do proprio cenario
addpath(fileparts(thisdir));  % raiz: SimRIS_v18
cd(thisdir);                  % le/grava resultados na pasta do cenario
% ------------------------------------------------------------------------

disp('=================================================================');
disp('OTIMIZACAO CENARIO B - 2 Tx ASSIMETRICOS + OBSTACULOS');
disp('=================================================================');

cfg = cenarioB_config();
pontos_rx = gera_usuarios_cenarioB(cfg);
fprintf('Usuarios sorteados: %d (semente %d)\n', size(pontos_rx,1), cfg.seed_users);

poolobj = gcp('nocreate');
if isempty(poolobj)
    disp('A iniciar Parallel Pool...');
    parpool;
end

N_lista = [64, 256];

out = 'convergencia_cenarioB';
if ~exist(out, 'dir'); mkdir(out); end

% guarda os usuarios sorteados para reuso (boxplot/topologia)
save(fullfile(out, 'usuarios_cenarioB.mat'), 'pontos_rx', 'cfg');

csv_path = fullfile(out, 'progresso_cenarioB.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, 'N,iteracao,fitness,x,y,z\n');
fclose(fid);
fprintf('CSV de progresso: %s\n', csv_path);

resultado = struct();
tic;
for kn = 1:numel(N_lista)
    N = N_lista(kn);
    fprintf('\n>>> OTIMIZANDO RIS (CENARIO B): N = %d <<<\n', N);

    [best, fbest, conv, sc] = fpa_cenarioB([20, 0.8], N, pontos_rx, cfg, csv_path);
    [~, taxas, taxas_base] = calculamediavazao_cenarioB(best, sc, N, pontos_rx, cfg);

    parede = 'Sul'; if best(2) >= cfg.Ub(2)/2; parede = 'Norte'; end
    fprintf('OK! Fitness: %.3f | RIS media: %.2f min: %.2f | Baseline media: %.2f min: %.2f\n', ...
            fbest, mean(taxas), min(taxas), mean(taxas_base), min(taxas_base));
    fprintf('    Melhor posicao: [%.1f, %.1f, %.1f] (parede %s)\n', ...
            best(1), best(2), best(3), parede);

    resultado(kn).N = N;
    resultado(kn).best = best;
    resultado(kn).fbest = fbest;
    resultado(kn).conv = conv;
    resultado(kn).parede = parede;
    resultado(kn).media = mean(taxas);
    resultado(kn).minimo = min(taxas);
    resultado(kn).media_base = mean(taxas_base);
    resultado(kn).minimo_base = min(taxas_base);

    fig = figure('Visible', 'off', 'Position', [100 100 760 520]);
    plot(1:numel(conv), conv, '-o', 'LineWidth', 1.6, 'MarkerSize', 3, ...
         'Color', [0.20 0.45 0.85], 'MarkerFaceColor', [0.20 0.45 0.85]);
    xlabel('Iteration', 'FontSize', 12);
    ylabel('Best Fitness (0.7\cdotmean + 0.3\cdotmin)', 'FontSize', 12);
    grid on; box on;
    fn = fullfile(out, sprintf('convergencia_cenarioB_N%d.png', N));
    saveas(fig, fn); close(fig);
    fprintf('    Gerado: %s\n', fn);
end
tempo = toc;

save(fullfile(out, 'resultado_cenarioB.mat'), 'resultado', 'cfg', 'pontos_rx');

disp('=================================================================');
fprintf('OTIMIZACAO CENARIO B CONCLUIDA EM %.2f SEGUNDOS.\n', tempo);
