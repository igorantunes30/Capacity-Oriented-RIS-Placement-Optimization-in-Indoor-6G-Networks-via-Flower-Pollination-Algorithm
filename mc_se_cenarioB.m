% =========================================================================
% mc_se_cenarioB.m
% 15 realizacoes Monte Carlo de SE por UE (96 UEs) p/ 3 cenarios:
%   - Baseline (sem RIS, com interferencia)
%   - RIS N=64  na melhor posicao do CSV de otimizacao
%   - RIS N=256 na melhor posicao do CSV de otimizacao
%
% NAO faz media MC dentro do calculo: salva uma linha por (realizacao, UE)
% no CSV. Modelo SINR multi-Tx (servidor + interferente via mesma config Phi).
% Roda SERIAL (sem parpool). Saida: mc_se_ues_cenarioB.csv
% =========================================================================
clear; clc; close all;
rng(1);

R_mc = 15;                               % numero de realizacoes
thisdir = fileparts(mfilename('fullpath'));
addpath(thisdir);
addpath(fileparts(thisdir));             % pasta cenarioB_obstaculos
addpath(fileparts(fileparts(thisdir)));  % raiz RecSurfaces (SimRIS_v18)
cd(thisdir);

cfg = cenarioB_config();
pos_Rx = gera_usuarios_cenarioB(cfg);
nU = size(pos_Rx, 1);

% --- melhor posicao por N do CSV de otimizacao --------------------------
csv_path = fullfile('..', 'convergencia_cenarioB', 'progresso_cenarioB.csv');
T = readmatrix(csv_path);                % N, iter, fitness, x, y, z
Ns = [64, 256];
best = containers.Map('KeyType','double','ValueType','any');
for N = Ns
    sub = T(T(:,1) == N, :);
    [fbest, ~] = max(sub(:,3));
    cands = sub(sub(:,3) == fbest, :);
    [~, k] = max(cands(:,2));
    best(N) = cands(k, 4:6);
    fprintf('N=%3d : best pos [%.2f, %.2f, %.2f] | fitness %.4f\n', ...
            N, cands(k,4), cands(k,5), cands(k,6), fbest);
end

Tx1 = cfg.Tx1; Tx2 = cfg.Tx2; L = cfg.Lsala;
SNR = cfg.SNR_lin; B = cfg.blockers;
Env = cfg.Env; Fr = cfg.Freq; AT = cfg.ArrayType; Nt = cfg.Nt; Nr = cfg.Nr;
sc = 1;                                  % paredes Sul/Norte -> Scenario 1
a_lin = @(dB) 10^(-dB/20);

pos64  = best(64);
pos256 = best(256);

% --- oclusao por UE (geometria fixa) ------------------------------------
% por UE: aD1, aD2 (baseline) + (aH1,aH2,aG_pos64) + (aH1_256,aH2_256,aG_pos256)
aD1 = zeros(nU,1); aD2 = zeros(nU,1);
aH1_64 = zeros(nU,1); aH2_64 = zeros(nU,1); aG_64 = zeros(nU,1);
aH1_256= zeros(nU,1); aH2_256= zeros(nU,1); aG_256= zeros(nU,1);
for i = 1:nU
    Rx = pos_Rx(i,:);
    aD1(i)    = a_lin(occlusao_att_dB(Tx1, Rx, B));
    aD2(i)    = a_lin(occlusao_att_dB(Tx2, Rx, B));
    aH1_64(i) = a_lin(occlusao_att_dB(Tx1, pos64,  B));
    aH2_64(i) = a_lin(occlusao_att_dB(Tx2, pos64,  B));
    aG_64(i)  = a_lin(occlusao_att_dB(pos64, Rx,   B));
    aH1_256(i)= a_lin(occlusao_att_dB(Tx1, pos256, B));
    aH2_256(i)= a_lin(occlusao_att_dB(Tx2, pos256, B));
    aG_256(i) = a_lin(occlusao_att_dB(pos256, Rx,  B));
end

% --- buffers de saida ---------------------------------------------------
SE_BL   = nan(R_mc, nU);
SE_64   = nan(R_mc, nU);
SE_256  = nan(R_mc, nU);

% --- frames espelhados (constantes p/ Tx2 em SimRIS) --------------------
Tx2_m   = [L - Tx2(1),   Tx2(2),   Tx2(3)];
pos64_m = [L - pos64(1), pos64(2), pos64(3)];
pos256_m= [L - pos256(1),pos256(2),pos256(3)];

t0 = tic;
for r = 1:R_mc
    for i = 1:nU
        Rx   = pos_Rx(i,:);
        Rx_m = [L - Rx(1), Rx(2), Rx(3)];
        try
            % ===== draw N=64 =====
            [H1_64, G1_64, D1]    = SimRIS_v18(Env, sc, Fr, AT, 64, Nt, Nr, Tx1,   Rx,   pos64);
            [H2_64, G2_64, D2]    = SimRIS_v18(Env, sc, Fr, AT, 64, Nt, Nr, Tx2_m, Rx_m, pos64_m);

            % baseline (so D, com oclusao) -- usa este draw p/ D1,D2
            D1bl = D1 * aD1(i); D2bl = D2 * aD2(i);
            potD1 = abs(D1bl)^2; potD2 = abs(D2bl)^2;
            potD_sig = max(potD1, potD2);
            potD_int = min(potD1, potD2);
            SINR_bl  = (SNR * potD_sig) / (SNR * potD_int + 1);
            SE_BL(r,i) = log2(1 + SINR_bl);

            % RIS N=64 em pos64
            Hv1 = H1_64(:)*aH1_64(i); Gv1 = G1_64(:)*aG_64(i); D1c = D1*aD1(i);
            Hv2 = H2_64(:)*aH2_64(i); Gv2 = G2_64(:)*aG_64(i); D2c = D2*aD2(i);
            th1 = angle(D1c) - (angle(Hv1) + angle(Gv1));
            th2 = angle(D2c) - (angle(Hv2) + angle(Gv2));
            C1  = (transpose(Gv1)*diag(exp(1j*th1))*Hv1) + D1c;
            C2  = (transpose(Gv2)*diag(exp(1j*th2))*Hv2) + D2c;
            pot1 = abs(C1)^2; pot2 = abs(C2)^2;
            if pot1 >= pot2
                Cint = (transpose(Gv2)*diag(exp(1j*th1))*Hv2) + D2c;
                pot_sig = pot1; pot_int = abs(Cint)^2;
            else
                Cint = (transpose(Gv1)*diag(exp(1j*th2))*Hv1) + D1c;
                pot_sig = pot2; pot_int = abs(Cint)^2;
            end
            SINR = (SNR*pot_sig) / (SNR*pot_int + 1);
            SE_64(r,i) = log2(1 + SINR);

            % ===== draw N=256 em pos256 =====
            [H1_256, G1_256, ~] = SimRIS_v18(Env, sc, Fr, AT, 256, Nt, Nr, Tx1,   Rx,   pos256);
            [H2_256, G2_256, ~] = SimRIS_v18(Env, sc, Fr, AT, 256, Nt, Nr, Tx2_m, Rx_m, pos256_m);

            Hv1 = H1_256(:)*aH1_256(i); Gv1 = G1_256(:)*aG_256(i);
            Hv2 = H2_256(:)*aH2_256(i); Gv2 = G2_256(:)*aG_256(i);
            % D1,D2 reaproveitados do draw N=64 (mesmas posicoes Tx-Rx)
            th1 = angle(D1c) - (angle(Hv1) + angle(Gv1));
            th2 = angle(D2c) - (angle(Hv2) + angle(Gv2));
            C1  = (transpose(Gv1)*diag(exp(1j*th1))*Hv1) + D1c;
            C2  = (transpose(Gv2)*diag(exp(1j*th2))*Hv2) + D2c;
            pot1 = abs(C1)^2; pot2 = abs(C2)^2;
            if pot1 >= pot2
                Cint = (transpose(Gv2)*diag(exp(1j*th1))*Hv2) + D2c;
                pot_sig = pot1; pot_int = abs(Cint)^2;
            else
                Cint = (transpose(Gv1)*diag(exp(1j*th2))*Hv1) + D1c;
                pot_sig = pot2; pot_int = abs(Cint)^2;
            end
            SINR = (SNR*pot_sig) / (SNR*pot_int + 1);
            SE_256(r,i) = log2(1 + SINR);

        catch ME
            warning('Erro r=%d UE=%d : %s', r, i, ME.message);
        end
    end
    fprintf('  realizacao %2d/%2d done | acumulado %.1f min\n', ...
            r, R_mc, toc(t0)/60);
end
fprintf('\nMC concluido em %.1f min (R=%d, UEs=%d)\n', toc(t0)/60, R_mc, nU);

% --- CSV long: realization, UE, baseline, N64, N256 ---------------------
fn_csv = 'mc_se_ues_cenarioB.csv';
fid = fopen(fn_csv, 'w');
fprintf(fid, 'realizacao,UE,baseline,RIS_N64,RIS_N256\n');
for r = 1:R_mc
    for i = 1:nU
        fprintf(fid, '%d,%d,%.6f,%.6f,%.6f\n', r, i, ...
                SE_BL(r,i), SE_64(r,i), SE_256(r,i));
    end
end
fclose(fid);
fprintf('Gerado: %s (%d linhas)\n', fn_csv, R_mc*nU);

% --- .mat p/ reuso -------------------------------------------------------
save('mc_se_ues_cenarioB.mat', 'SE_BL', 'SE_64', 'SE_256', ...
     'pos64', 'pos256', 'R_mc', 'cfg', 'pos_Rx');
fprintf('Gerado: mc_se_ues_cenarioB.mat\n');

% --- resumo --------------------------------------------------------------
fprintf('\nResumo SE (bps/Hz) - %d amostras por grupo (R*UE):\n', R_mc*nU);
fprintf('%-15s %8s %8s %8s %8s %8s %8s\n','Grupo','min','P5','P25','mediana','media','max');
print_stat = @(name, V) fprintf('%-15s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', ...
    name, min(V(:)), quantile_local(V(:),0.05), quantile_local(V(:),0.25), ...
    quantile_local(V(:),0.50), mean(V(:),'omitnan'), max(V(:)));
print_stat('Baseline',  SE_BL);
print_stat('RIS N=64',  SE_64);
print_stat('RIS N=256', SE_256);

function q = quantile_local(v, p)
    v = sort(v(~isnan(v)));
    n = numel(v); if n==0; q=NaN; return; end
    pos = p*n + 0.5; pos = min(max(pos,1),n);
    lo = floor(pos); hi = ceil(pos);
    if lo==hi; q=v(lo); else; q = v(lo)+(pos-lo)*(v(hi)-v(lo)); end
end
