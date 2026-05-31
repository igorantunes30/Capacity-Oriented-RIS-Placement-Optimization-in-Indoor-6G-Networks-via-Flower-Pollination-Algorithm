% =========================================================================
% teste_ris_alta_cenarioB.m
% Igual ao teste_ris_fixa, mas com a RIS ALTA em [x=40, Sul (y=0), z=2.9],
% acima do topo dos obstaculos (2.4 m) -> caminho Tx->RIS LIVRE. Contraste
% direto com a RIS baixa (z=1.5): mesma posicao lateral, so muda a altura.
% Roda SERIAL (sem parpool) -> nao interfere na otimizacao em andamento.
% =========================================================================
clear; clc; close all;
rng(1);

% --- bootstrap de caminhos (cenario em pasta propria) -------------------
thisdir = fileparts(mfilename('fullpath'));
addpath(thisdir);
addpath(fileparts(thisdir));  % raiz: SimRIS_v18
cd(thisdir);
% ------------------------------------------------------------------------

cfg = cenarioB_config();
pos_Rx = gera_usuarios_cenarioB(cfg);

pos_RIS = [40, 0, 2.9];   % x=40 na parede Sul (y=0), altura ALTA z=2.9
sc = 1;                   % parede lateral (Sul) -> Scenario 1
N_list = [64, 256];

fprintf('=================================================================\n');
fprintf('TESTE FIXO CENARIO B - RIS ALTA (sem otimizacao) | %d usuarios | R_mc=%d\n', ...
        size(pos_Rx,1), cfg.R_mc);
fprintf('RIS fixa em [%.0f, %.0f, %.1f] (parede Sul, acima dos obstaculos)\n', pos_RIS);
fprintf('=================================================================\n');

grupos = {}; rotulos = {}; resumo = struct();
base_feito = false; kk = 0;

for N = N_list
    fprintf('\n>>> Avaliando N = %d ...\n', N);
    [~, taxas, taxas_base] = calculamediavazao_cenarioB(pos_RIS, sc, N, pos_Rx, cfg);

    if ~base_feito
        grupos{end+1} = taxas_base(:); rotulos{end+1} = 'Baseline (sem RIS)'; %#ok<SAGROW>
        fprintf('  SEM RIS  : mediana %.3f | media %.3f | min %.3f\n', ...
                median(taxas_base), mean(taxas_base), min(taxas_base));
        resumo.base = [median(taxas_base), mean(taxas_base), min(taxas_base)];
        base_feito = true;
    end

    grupos{end+1} = taxas(:); rotulos{end+1} = sprintf('RIS N=%d', N); %#ok<SAGROW>
    fprintf('  COM RIS  : mediana %.3f | media %.3f | min %.3f\n', ...
            median(taxas), mean(taxas), min(taxas));
    kk = kk + 1;
    resumo.ris(kk,:) = [N, median(taxas), mean(taxas), min(taxas)];
end

save('teste_ris_alta_cenarioB.mat', 'pos_RIS', 'sc', 'N_list', 'grupos', 'rotulos', 'resumo', 'cfg');

% --- Boxplot comparativo (manual, sem Statistics Toolbox) ----------------
fig = figure('Visible','off','Position',[100 100 820 560]); hold on;
cores = [0.70 0.70 0.75; 0.20 0.45 0.85; 0.85 0.30 0.30];
larg = 0.5; ymax = 0;
for g = 1:numel(grupos)
    v = sort(grupos{g});
    q1 = pct(v,25); med = pct(v,50); q3 = pct(v,75); iqr = q3-q1;
    lo = max(min(v), q1-1.5*iqr); hi = min(max(v), q3+1.5*iqr);
    outl = v(v<lo | v>hi);
    cor = cores(mod(g-1,size(cores,1))+1,:);
    xL = g-larg/2; xR = g+larg/2;
    fill([xL xR xR xL],[q1 q1 q3 q3], cor, 'FaceAlpha',0.6,'EdgeColor','k','LineWidth',1);
    plot([xL xR],[med med],'k-','LineWidth',2);
    hM = plot(g, mean(grupos{g}),'d','MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','k','MarkerSize',7);
    plot([g g],[q3 hi],'k--','LineWidth',1); plot([g g],[lo q1],'k--','LineWidth',1);
    plot([g-0.15 g+0.15],[hi hi],'k-'); plot([g-0.15 g+0.15],[lo lo],'k-');
    if ~isempty(outl); plot(g*ones(size(outl)), outl,'r+','MarkerSize',5); end
    ymax = max(ymax, max(v));
end
hold off;
set(gca,'XTick',1:numel(rotulos),'XTickLabel',rotulos,'FontSize',11);
xlim([0.4 numel(rotulos)+0.6]); ylim([0 max(ymax*1.15, 0.1)]);
ylabel('Eficiencia Espectral (bps/Hz)','FontSize',12);
title({'Cenario B - Teste com RIS ALTA (sem otimizacao)', ...
       sprintf('RIS em [%.0f, %.0f, %.1f] (Sul, acima dos obstaculos) | %d dB | %d usuarios', ...
               pos_RIS(1),pos_RIS(2),pos_RIS(3), cfg.att_obstaculo, size(pos_Rx,1))}, ...
      'FontSize',12);
legend(hM, {'Media'}, 'Location','northeast','FontSize',10);
grid on; box on;

fn = 'teste_ris_alta_cenarioB.png';
saveas(fig, fn); close(fig);
fprintf('\nGerado: %s\n', fn);

function q = pct(v, p)
    n = numel(v);
    if n==1; q=v(1); return; end
    pos = (p/100)*n + 0.5; pos = min(max(pos,1),n);
    lo = floor(pos); hi = ceil(pos);
    if lo==hi; q=v(lo); else; q = v(lo)+(pos-lo)*(v(hi)-v(lo)); end
end
