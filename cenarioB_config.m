function cfg = cenarioB_config()
% Configuracao do CENARIO B (sala assimetrica, paredes grossas, obstaculos
% bloqueando a visada direta dos Tx). Tudo em coordenadas REAIS da sala.
% Difere do cenario 2 Tx original por:
%   - Tx em alturas/y assimetricos (quebra a simetria em x);
%   - obstaculos (caixas) na frente de cada Tx => bloqueiam LOS direto;
%   - uma parede interna GROSSA dividindo o salao;
%   - usuarios sorteados aleatoriamente (semente fixa).

    cfg = struct();

    % --- Sala (mantida 75x50x3.5 p/ coerencia com os bounds do SimRIS) ---
    cfg.room = [75, 50, 3.5];

    % --- Transmissores ASSIMETRICOS (y diferentes => sem simetria em x) ---
    cfg.Tx1 = [0,  15, 2.5];     % parede Oeste, mais ao Sul
    cfg.Tx2 = [75, 40, 2.5];     % parede Leste, mais ao Norte
    cfg.Lsala = cfg.room(1);     % usado no referencial espelhado do Tx2

    % --- Radio / SimRIS ---
    cfg.Pt_dBm = 30; cfg.Pn_dBm = -80;
    cfg.SNR_lin = 10^((cfg.Pt_dBm - cfg.Pn_dBm)/10);
    cfg.Env = 1; cfg.Freq = 28; cfg.ArrayType = 2; cfg.Nt = 1; cfg.Nr = 1;
    cfg.R_mc = 20;
    cfg.alpha = 0.7;             % fitness = alpha*media + (1-alpha)*min

    % --- Bloqueadores: [xmin xmax ymin ymax zmin zmax att_dB] ---
    % Apenas "telas" (obstaculos) na frente de cada Tx -- SEM parede interna.
    % Topo dos obstaculos em 2.4 m (abaixo do Tx em 2.5 m e da RIS >=2.6 m):
    % o direto Tx->Rx (desce ate 1 m) atravessa e e atenuado, mas o caminho
    % Tx->RIS passa POR CIMA (>=2.5 m) e fica livre.
    cfg.att_obstaculo = 30;      % atenuacao do obstaculo (dB)
    cfg.z_obst = 2.4;            % topo dos obstaculos
    cfg.blockers = [
        % tela na frente do Tx1 (assimetrica: banda Sul), topo 2.4 m
       16.0, 19.0,  4.0, 30.0, 0.0, 2.4, cfg.att_obstaculo;
        % tela na frente do Tx2 (assimetrica: banda Norte), topo 2.4 m
       56.0, 59.0, 20.0, 46.0, 0.0, 2.4, cfg.att_obstaculo;
    ];

    % --- Usuarios aleatorios ---
    cfg.num_users = 96;
    cfg.seed_users = 7;          % semente fixa p/ reprodutibilidade
    cfg.user_z = 1.0;            % altura do receptor
    cfg.margem = 2.0;            % afastamento minimo das paredes
    cfg.dist_min_tx = 3.0;       % afastamento minimo dos Tx

    % --- Espaco de busca da RIS (paredes laterais Sul/Norte) ---
    % z >= 2.6 m: RIS sempre acima do topo dos obstaculos (2.4) e do Tx (2.5),
    % garantindo que o caminho Tx->RIS passe por cima das telas.
    cfg.Lb = [0, 0, 2.6];
    cfg.Ub = [75, 50, 3.0];
end
