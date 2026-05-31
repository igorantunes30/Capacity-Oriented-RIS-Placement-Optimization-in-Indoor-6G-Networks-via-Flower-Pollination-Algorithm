% =========================================================================
% snr_por_ue_cenarioB.m
% SNR (sem interferencia) por UE usando as POSICOES OTIMAS da RIS (N=64 e
% N=256) lidas de ../convergencia_cenarioB/progresso_cenarioB.csv.
%
% Diferenca p/ calculamediavazao_cenarioB (que usa SINR multi-Tx):
%   aqui o "servidor" e o Tx com maior |C|^2 cascateado, mas o outro Tx NAO
%   entra como interferente. So ruido:
%       SNR = SNR_lin * |C_sig|^2          (linear)
%       SNR_dB = 10*log10(SNR)
%   Baseline (sem RIS): servidor = direto mais forte, SNR = SNR_lin*max(|D|^2).
%
% Media Monte Carlo: media de |C_sig|^2 (potencia) sobre R realizacoes, depois
% converte p/ SNR. Roda SERIAL (sem parpool). Saida: snr_por_ue_cenarioB.csv
% =========================================================================
clear; clc; close all;
rng(1);

thisdir = fileparts(mfilename('fullpath'));
addpath(thisdir);
addpath(fileparts(thisdir));             % cenarioB_obstaculos
addpath(fileparts(fileparts(thisdir)));  % raiz RecSurfaces (SimRIS_v18)
cd(thisdir);

cfg = cenarioB_config();
pos_Rx = gera_usuarios_cenarioB(cfg);
nU = size(pos_Rx, 1);

Tx1 = cfg.Tx1; Tx2 = cfg.Tx2; L = cfg.Lsala;
SNR = cfg.SNR_lin; B = cfg.blockers;
Env = cfg.Env; Fr = cfg.Freq; AT = cfg.ArrayType; Nt = cfg.Nt; Nr = cfg.Nr;
R = cfg.R_mc;
sc = 1;                                  % paredes Sul/Norte -> Scenario 1
a_lin = @(dB) 10^(-dB/20);

% --- melhor posicao por N do CSV de otimizacao --------------------------
csv_path = fullfile('..', 'convergencia_cenarioB', 'progresso_cenarioB.csv');
T = readmatrix(csv_path);                % N, iter, fitness, x, y, z
Ns = [64, 256];
best = containers.Map('KeyType','double','ValueType','any');
for N = Ns
    sub = T(T(:,1) == N, :);
    [fbest, ~] = max(sub(:,3));
    cands = sub(sub(:,3) == fbest, :);
    [~, k] = max(cands(:,2));            % em empate, ultima iteracao
    best(N) = cands(k, 4:6);
    fprintf('N=%3d : best pos [%.2f, %.2f, %.2f] | fitness %.4f\n', ...
            N, cands(k,4), cands(k,5), cands(k,6), fbest);
end

% --- buffers ------------------------------------------------------------
snr_bl  = zeros(nU,1);                    % baseline (independe de N)
snr_64  = zeros(nU,1);
snr_256 = zeros(nU,1);

t0 = tic;
for kn = 1:numel(Ns)
    N = Ns(kn);
    pos_RIS = best(N);
    Tx2_m   = [L - Tx2(1),       Tx2(2),       Tx2(3)];
    RIS_m   = [L - pos_RIS(1),   pos_RIS(2),   pos_RIS(3)];

    for i = 1:nU
        Rx   = pos_Rx(i,:);
        Rx_m = [L - Rx(1), Rx(2), Rx(3)];

        aD1 = a_lin(occlusao_att_dB(Tx1, Rx,      B));
        aH1 = a_lin(occlusao_att_dB(Tx1, pos_RIS, B));
        aG  = a_lin(occlusao_att_dB(pos_RIS, Rx,  B));
        aD2 = a_lin(occlusao_att_dB(Tx2, Rx,      B));
        aH2 = a_lin(occlusao_att_dB(Tx2, pos_RIS, B));

        acc_sig = 0; acc_Dsig = 0; ok = 0;
        for r = 1:R
            try
                [H1, G1, D1] = SimRIS_v18(Env, sc, Fr, AT, N, Nt, Nr, Tx1,   Rx,   pos_RIS);
                Hv1 = H1(:)*aH1; Gv1 = G1(:)*aG; D1 = D1*aD1;
                th1 = angle(D1) - (angle(Hv1) + angle(Gv1));
                C1  = (transpose(Gv1)*diag(exp(1j*th1))*Hv1) + D1;
                pot1 = abs(C1)^2; potD1 = abs(D1)^2;

                [H2, G2, D2] = SimRIS_v18(Env, sc, Fr, AT, N, Nt, Nr, Tx2_m, Rx_m, RIS_m);
                Hv2 = H2(:)*aH2; Gv2 = G2(:)*aG; D2 = D2*aD2;
                th2 = angle(D2) - (angle(Hv2) + angle(Gv2));
                C2  = (transpose(Gv2)*diag(exp(1j*th2))*Hv2) + D2;
                pot2 = abs(C2)^2; potD2 = abs(D2)^2;

                acc_sig  = acc_sig  + max(pot1, pot2);     % servidor c/ RIS
                acc_Dsig = acc_Dsig + max(potD1, potD2);   % servidor baseline
                ok = ok + 1;
            catch ME
                warning('N=%d UE=%d r=%d : %s', N, i, r, ME.message);
            end
        end
        if ok > 0
            snr_ris = SNR * (acc_sig / ok);
            if N == 64;  snr_64(i)  = snr_ris; end
            if N == 256; snr_256(i) = snr_ris; end
            if kn == 1;  snr_bl(i)  = SNR * (acc_Dsig / ok); end  % baseline 1x
        end
    end
    fprintf('  N=%d done | %.1f min\n', N, toc(t0)/60);
end
fprintf('\nConcluido em %.1f min (R=%d, UEs=%d)\n', toc(t0)/60, R, nU);

% --- CSV: UE, SNR linear + dB por grupo ---------------------------------
to_dB = @(x) 10*log10(x);
fn = fullfile('..', 'convergencia_cenarioB', 'snr_por_ue_cenarioB.csv');
fid = fopen(fn, 'w');
fprintf(fid, ['UE,snr_baseline_lin,snr_baseline_dB,' ...
              'snr_N64_lin,snr_N64_dB,snr_N256_lin,snr_N256_dB\n']);
for i = 1:nU
    fprintf(fid, '%d,%.6e,%.4f,%.6e,%.4f,%.6e,%.4f\n', i, ...
        snr_bl(i),  to_dB(snr_bl(i)), ...
        snr_64(i),  to_dB(snr_64(i)), ...
        snr_256(i), to_dB(snr_256(i)));
end
fclose(fid);
fprintf('Gerado: %s (%d UEs)\n', fn, nU);

% --- resumo (dB) --------------------------------------------------------
fprintf('\nResumo SNR (dB):\n');
fprintf('%-15s %8s %8s %8s %8s\n','Grupo','min','mediana','media','max');
prn = @(name,v) fprintf('%-15s %8.2f %8.2f %8.2f %8.2f\n', name, ...
    min(to_dB(v)), median(to_dB(v)), mean(to_dB(v)), max(to_dB(v)));
prn('Baseline',  snr_bl);
prn('RIS N=64',  snr_64);
prn('RIS N=256', snr_256);
