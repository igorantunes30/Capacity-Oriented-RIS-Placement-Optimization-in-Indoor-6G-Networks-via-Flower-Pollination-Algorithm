function [fitness, taxas, taxas_baseline] = calculamediavazao_cenarioB(u, cenario_simris, N_elementos, pos_Rx, cfg)
% Vazao no CENARIO B (2 Tx assimetricos + obstaculos + parede grossa).
% Modelo SINR multi-Tx: o Tx mais forte (com RIS) e o SERVIDOR; o outro Tx
% e INTERFERENTE e contribui com seu canal cascateado avaliado sob a MESMA
% configuracao de fase da RIS (theta otimizado p/ o servidor) + direto.
% A oclusao geometrica adiciona perda (dB) aos caminhos que cruzam
% obstaculos/paredes: direto (D), Tx->RIS (H) e RIS->Rx (G). O Tx2 (parede
% Leste) usa referencial ESPELHADO em x apenas para CHAMAR o SimRIS; a
% oclusao e sempre calculada em coordenadas REAIS.
%
% Eficiencia espectral por usuario:
%   SE = log2(1 + SINR),  SINR = (SNR_lin * |C_sig|^2) /
%                                (SNR_lin * |C_int|^2 + 1)
% onde SNR_lin = Pt/Pn (linear). O "+1" no denominador representa a
% potencia de ruido normalizada (signal e interferente sao escalados pelo
% mesmo Pt, e o ruido pelo mesmo Pn => razao SINR independente das
% unidades absolutas).
%
% Saidas:
%   fitness        = alpha*mean(taxas) + (1-alpha)*min(taxas)   [com RIS]
%   taxas          = SE por usuario COM a RIS (bps/Hz)
%   taxas_baseline = SE por usuario SEM a RIS (so diretos, com interf.)

    if nargin < 5; cfg = cenarioB_config(); end

    Tx1 = cfg.Tx1; Tx2 = cfg.Tx2; L = cfg.Lsala;
    SNR = cfg.SNR_lin; B = cfg.blockers;
    Env = cfg.Env; Fr = cfg.Freq; AT = cfg.ArrayType; Nt = cfg.Nt; Nr = cfg.Nr;
    R = cfg.R_mc;

    num = size(pos_Rx, 1);
    taxas = zeros(num, 1);
    taxas_baseline = zeros(num, 1);
    pos_RIS = [u(1), u(2), u(3)];

    % atenuacoes lineares (amplitude) dependentes so da geometria (fixas no MC)
    a_lin = @(dB) 10^(-dB/20);

    for i = 1:num
        Rx = pos_Rx(i,:);

        % --- oclusao (coordenadas reais) ---
        aD1 = a_lin(occlusao_att_dB(Tx1, Rx,      B));   % direto Tx1
        aH1 = a_lin(occlusao_att_dB(Tx1, pos_RIS, B));   % Tx1 -> RIS
        aG1 = a_lin(occlusao_att_dB(pos_RIS, Rx,  B));   % RIS  -> Rx
        aD2 = a_lin(occlusao_att_dB(Tx2, Rx,      B));   % direto Tx2
        aH2 = a_lin(occlusao_att_dB(Tx2, pos_RIS, B));   % Tx2 -> RIS
        % G (RIS->Rx) e comum aos dois Tx (mesma RIS, mesmo Rx)
        aG  = aG1;

        acc = 0; accB = 0; ok = 0;
        for r = 1:R
            try
                % --- Tx1 -> (RIS alinhada a Tx1) -> Rx ---
                [H1, G1, D1] = SimRIS_v18(Env, cenario_simris, Fr, AT, ...
                    N_elementos, Nt, Nr, Tx1, Rx, pos_RIS);
                Hv1 = H1(:) * aH1;  Gv1 = G1(:) * aG;  D1 = D1 * aD1;
                th1 = angle(D1) - (angle(Hv1) + angle(Gv1));
                C1  = (transpose(Gv1) * diag(exp(1j*th1)) * Hv1) + D1;
                pot1 = abs(C1)^2;  potD1 = abs(D1)^2;

                % --- Tx2 -> (RIS alinhada a Tx2) -> Rx  (frame espelhado) ---
                Tx2_m = [L - Tx2(1), Tx2(2), Tx2(3)];
                RIS_m = [L - pos_RIS(1), pos_RIS(2), pos_RIS(3)];
                Rx_m  = [L - Rx(1), Rx(2), Rx(3)];
                [H2, G2, D2] = SimRIS_v18(Env, cenario_simris, Fr, AT, ...
                    N_elementos, Nt, Nr, Tx2_m, Rx_m, RIS_m);
                Hv2 = H2(:) * aH2;  Gv2 = G2(:) * aG;  D2 = D2 * aD2;
                th2 = angle(D2) - (angle(Hv2) + angle(Gv2));
                C2  = (transpose(Gv2) * diag(exp(1j*th2)) * Hv2) + D2;
                pot2 = abs(C2)^2;  potD2 = abs(D2)^2;

                % RIS alinha ao Tx mais forte (servidor). O outro Tx vira
                % INTERFERENTE: seu canal cascateado e reavaliado sob a
                % MESMA configuracao de fase da RIS (theta do servidor) e
                % somado ao seu direto.
                if pot1 >= pot2
                    % Servidor = Tx1 (theta = th1). Interferente = Tx2.
                    Cint    = (transpose(Gv2) * diag(exp(1j*th1)) * Hv2) + D2;
                    pot_sig = pot1;
                    pot_int = abs(Cint)^2;
                else
                    % Servidor = Tx2 (theta = th2). Interferente = Tx1.
                    Cint    = (transpose(Gv1) * diag(exp(1j*th2)) * Hv1) + D1;
                    pot_sig = pot2;
                    pot_int = abs(Cint)^2;
                end

                % SINR normalizado (Pt cancela no numerador/denominador):
                %   SINR = (SNR_lin * |C_sig|^2) / (SNR_lin * |C_int|^2 + 1)
                SINR    = (SNR * pot_sig) / (SNR * pot_int + 1);

                % Baseline sem RIS: servidor = Tx do direto mais forte,
                % interferente = direto do outro Tx.
                potD_sig = max(potD1, potD2);
                potD_int = min(potD1, potD2);
                SINR_bl  = (SNR * potD_sig) / (SNR * potD_int + 1);

                acc  = acc  + log2(1 + SINR);
                accB = accB + log2(1 + SINR_bl);
                ok = ok + 1;
            catch ME
                disp(['Erro SimRIS/Algebra (cenarioB): ' ME.message]);
            end
        end
        if ok > 0
            taxas(i)          = acc  / ok;
            taxas_baseline(i) = accB / ok;
        end
    end

    fitness = cfg.alpha * mean(taxas) + (1 - cfg.alpha) * min(taxas);
end
