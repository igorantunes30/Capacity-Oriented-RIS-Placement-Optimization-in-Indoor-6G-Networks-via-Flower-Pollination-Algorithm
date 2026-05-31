% =========================================================================
% sanity_cenarioB.m
% Teste rapido do CENARIO B (sem otimizacao): valida o teste de oclusao e
% confirma que a RIS gera ganho sobre o baseline nos usuarios bloqueados,
% sem NaNs. Monte Carlo curto e poucos usuarios -> roda em segundos.
% =========================================================================
clear; clc;
rng(1);

% --- bootstrap de caminhos (cenario em pasta propria) -------------------
thisdir = fileparts(mfilename('fullpath'));
addpath(thisdir);             % funcoes do proprio cenario
addpath(fileparts(thisdir));  % raiz: SimRIS_v18
cd(thisdir);
% ------------------------------------------------------------------------

cfg = cenarioB_config();
cfg.R_mc = 6;     % Monte Carlo curto so para o sanity

% ---- 1) Teste direto da funcao de oclusao ----
disp('--- Teste de oclusao (att em dB) ---');
casos = {
    'Tx1(2.5) -> Rx atras da tela',  cfg.Tx1, [30, 15, 1.0];   % direto baixo: bloqueado
    'Tx1 -> Rx Norte (livre)',       cfg.Tx1, [30, 45, 1.0];   % fora da tela: livre
    'Tx1 -> RIS alta [37,50,2.9]',   cfg.Tx1, [37, 50, 2.9];   % passa por cima: livre
    'Tx2 -> Rx atras da tela Tx2',   cfg.Tx2, [50, 35, 1.0];   % direto baixo: bloqueado
};
for c = 1:size(casos,1)
    att = occlusao_att_dB(casos{c,2}, casos{c,3}, cfg.blockers);
    fprintf('  %-32s att = %2.0f dB\n', casos{c,1}, att);
end

% ---- 2) Usuarios de teste: alguns bloqueados, alguns livres ----
pos_teste = [
    30, 15, 1.0;   % atras da tela do Tx1 (direto bloqueado)
    30, 45, 1.0;   % Norte, mais livre
    65, 30, 1.0;   % atras da tela do Tx2
    50, 10, 1.0;   % Sul central
    50, 40, 1.0;   % Norte central
    50, 25, 1.0;   % centro
];

% RIS de teste alta na parede Norte (acima dos obstaculos de 2.6 m)
ris_teste = [37, 50, 2.9];
N = 64; sc = 1;

fprintf('\n--- Vazao por usuario: Baseline (sem RIS) vs RIS N=%d em [%.0f,%.0f,%.1f] ---\n', ...
        N, ris_teste(1), ris_teste(2), ris_teste(3));
[fit, taxas, taxas_base] = calculamediavazao_cenarioB(ris_teste, sc, N, pos_teste, cfg);
for i = 1:size(pos_teste,1)
    g = taxas(i) - taxas_base(i);
    fprintf('  Rx [%4.0f,%4.0f] | base %5.2f | RIS %5.2f | ganho %+5.2f | NaN=%d\n', ...
            pos_teste(i,1), pos_teste(i,2), taxas_base(i), taxas(i), g, ...
            isnan(taxas(i)) || isnan(taxas_base(i)));
end
fprintf('\nResumo: fitness(RIS)=%.3f | media base=%.2f | media RIS=%.2f | NaNs=%d\n', ...
        fit, mean(taxas_base), mean(taxas), sum(isnan(taxas)));
disp('Sanity-test concluido.');
