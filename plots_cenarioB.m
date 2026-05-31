% =========================================================================
% plots_cenarioB.m
% UNIFICA os 10 scripts de plotagem do Cenario B em um unico arquivo.
% Cada figura virou uma funcao local. Os helpers pct() e plot_box() sao
% compartilhados (antes duplicados em varios scripts).
%
% USO:
%   plots_cenarioB                      % roda TODAS as figuras
%   plots_cenarioB('topologia','cdf_sinr')   % roda so as escolhidas
%
% Nomes validos:
%   topologia            (era plot_topologia_cenarioB.m)
%   topologia_ris        (era plot_topologia_cenarioB_ris.m)
%   topologia_ris_otim   (era plot_topologia_cenarioB_ris_otim.m)
%   convergencia         (era plot_convergencia_cenarioB.m)
%   convergencia_N256    (era convergencia_N256.m)
%   boxplot_otim_parcial (era boxplot_otim_parcial.m)
%   boxplot_sinr         (era boxplot_sinr_cenarioB.m)
%   boxplot_mc           (era boxplot_mc_cenarioB.m)
%   cdf_mc               (era cdf_mc_cenarioB.m)
%   cdf_sinr             (era cdf_sinr_cenarioB.m)
%
% NOTA DE CAMINHOS: os caminhos de dados (CSV/MAT) foram mantidos como nos
% originais. Como os scripts viviam em 3 pastas diferentes (raiz,
% ris-interferencia, convergencia_cenarioB), os caminhos relativos assumem
% essas pastas. Ajuste cd/csv dentro de cada funcao se os dados estiverem
% em outro lugar.
% =========================================================================
function plots_cenarioB(varargin)
    all_figs = {'topologia','topologia_ris','topologia_ris_otim', ...
                'convergencia','convergencia_N256','boxplot_otim_parcial', ...
                'boxplot_sinr','boxplot_mc','cdf_mc','cdf_sinr'};
    if nargin == 0
        sel = all_figs;
    else
        sel = varargin;
    end
    for i = 1:numel(sel)
        switch sel{i}
            case 'topologia';            fig_topologia();
            case 'topologia_ris';        fig_topologia_ris();
            case 'topologia_ris_otim';   fig_topologia_ris_otim();
            case 'convergencia';         fig_convergencia();
            case 'convergencia_N256';    fig_convergencia_N256();
            case 'boxplot_otim_parcial'; fig_boxplot_otim_parcial();
            case 'boxplot_sinr';         fig_boxplot_sinr();
            case 'boxplot_mc';           fig_boxplot_mc();
            case 'cdf_mc';               fig_cdf_mc();
            case 'cdf_sinr';             fig_cdf_sinr();
            otherwise
                warning('Figura desconhecida: %s (ignorada)', sel{i});
        end
    end
end

% =========================================================================
% topologia (sem RIS) -- era plot_topologia_cenarioB.m
% =========================================================================
function fig_topologia()
    thisdir = fileparts(mfilename('fullpath'));
    addpath(thisdir);             % funcoes do proprio cenario
    addpath(fileparts(thisdir));  % raiz: SimRIS_v18
    cd(thisdir);                  % le/grava resultados na pasta do cenario

    cfg = cenarioB_config();
    pontos_rx = gera_usuarios_cenarioB(cfg);

    out = 'convergencia_cenarioB';
    if ~exist(out, 'dir'); mkdir(out); end
    save(fullfile(out, 'usuarios_cenarioB.mat'), 'pontos_rx', 'cfg');

    LX = cfg.room(1); LY = cfg.room(2);

    fig = figure('Visible','off','Position',[100 100 1000 680]);
    hold on;

    % --- contorno da sala (paredes Norte/Sul em verde, Leste/Oeste em preto) ---
    verde = [0 0.6 0];
    plot([0 0],   [0 LY], 'k-', 'LineWidth', 2);   % Oeste (x=0)
    plot([LX LX], [0 LY], 'k-', 'LineWidth', 2);   % Leste (x=LX)
    plot([0 LX],  [0 0],  '-', 'Color', verde, 'LineWidth', 2.5);   % Sul (y=0)
    plot([0 LX],  [LY LY],'-', 'Color', verde, 'LineWidth', 2.5);   % Norte (y=LY)

    % --- bloqueadores (telas-obstaculo em laranja) ---
    B = cfg.blockers;
    for k = 1:size(B,1)
        x1=B(k,1); x2=B(k,2); y1=B(k,3); y2=B(k,4);
        cor = [0.95 0.55 0.15];
        patch([x1 x2 x2 x1],[y1 y1 y2 y2], cor, 'FaceAlpha',0.55, ...
              'EdgeColor','k','LineWidth',1.2);
    end

    % --- usuarios ---
    hU = scatter(pontos_rx(:,1), pontos_rx(:,2), 26, [0.20 0.45 0.85], 'filled', ...
                 'MarkerEdgeColor','k','LineWidth',0.3);

    % --- transmissores ---
    hT = plot(cfg.Tx1(1), cfg.Tx1(2), 'p', 'MarkerSize',22, ...
              'MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor','k','LineWidth',1.2);
    plot(cfg.Tx2(1), cfg.Tx2(2), 'p', 'MarkerSize',22, ...
         'MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor','k','LineWidth',1.2);

    % --- posicoes otimas das RIS (N=64 e N=256) lidas da CSV de progresso -----
    Ns        = [64, 256];
    cores_ris = {[0.20 0.45 0.85], [0.85 0.30 0.30]};   % N=64 azul, N=256 vermelho
    ref_pos   = [59.2 50.0 3.0; 56.2 50.0 3.0];         % fallback (README)
    csv_path  = fullfile(out, 'progresso_cenarioB.csv');
    T = [];
    if exist(csv_path, 'file'); T = readmatrix(csv_path); end

    ris_handles = []; ris_labels = {};
    for kn = 1:numel(Ns)
        N = Ns(kn);
        pos = ref_pos(kn,:);
        if ~isempty(T)
            sub = T(T(:,1)==N, :);
            if ~isempty(sub)
                [~, ord] = sort(sub(:,2)); sub = sub(ord, :);
                pos = sub(end, 4:6);
            end
        end
        h = plot(pos(1), pos(2), 's', 'MarkerSize',16, ...
                 'MarkerFaceColor', cores_ris{kn}, 'MarkerEdgeColor','k','LineWidth',1.5);
        ris_handles(end+1) = h; %#ok<AGROW>
        ris_labels{end+1}  = sprintf('Optimal RIS N=%d', N); %#ok<AGROW>
    end

    % --- rotulos das paredes ---
    text(LX/2, -3, 'South Wall', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);
    text(LX/2, LY+3, 'North Wall', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);

    hold off; axis equal;
    xlim([-10 LX+12]); ylim([-8 LY+8]);
    xlabel('x (m)','FontSize',12); ylabel('y (m)','FontSize',12);
    legend([hT hU ris_handles], [{'Transmitters','Users'} ris_labels], ...
           'Location','northeastoutside','FontSize',10);
    grid on; box on;

    fn = 'topologia_cenarioB.png';
    exportgraphics(fig, fn, 'Resolution', 200); close(fig);   % corta margem branca
    fprintf('Gerado: %s | usuarios: %d\n', fn, size(pontos_rx,1));
end

% =========================================================================
% topologia_ris (RIS de teste, parede Sul) -- era plot_topologia_cenarioB_ris.m
% =========================================================================
function fig_topologia_ris()
    thisdir = fileparts(mfilename('fullpath'));
    addpath(thisdir);
    addpath(fileparts(thisdir));  % raiz: SimRIS_v18
    cd(thisdir);

    cfg = cenarioB_config();
    pontos_rx = gera_usuarios_cenarioB(cfg);

    ris_xy = [40, 0];            % x=40 na parede Sul (y=0)
    ris_zs = [1.5, 2.9];         % alturas testadas

    LX = cfg.room(1); LY = cfg.room(2);

    fig = figure('Visible','off','Position',[100 100 1000 680]);
    hold on;
    plot([0 LX LX 0 0], [0 0 LY LY 0], 'k-', 'LineWidth', 2);

    % --- obstaculos ---
    B = cfg.blockers;
    for k = 1:size(B,1)
        x1=B(k,1); x2=B(k,2); y1=B(k,3); y2=B(k,4); ztop=B(k,6); att=B(k,7);
        patch([x1 x2 x2 x1],[y1 y1 y2 y2], [0.95 0.55 0.15], 'FaceAlpha',0.55, ...
              'EdgeColor','k','LineWidth',1.2);
        text((x1+x2)/2, (y1+y2)/2, sprintf('Obstaculo\n%d dB (h=%.1fm)', att, ztop), ...
             'HorizontalAlignment','center','FontSize',9,'FontWeight','bold','Color',[0.1 0.1 0.1]);
    end

    % --- usuarios ---
    hU = scatter(pontos_rx(:,1), pontos_rx(:,2), 26, [0.20 0.45 0.85], 'filled', ...
                 'MarkerEdgeColor','k','LineWidth',0.3);

    % --- transmissores ---
    hT = plot(cfg.Tx1(1), cfg.Tx1(2), 'p', 'MarkerSize',22, ...
              'MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor','k','LineWidth',1.2);
    plot(cfg.Tx2(1), cfg.Tx2(2), 'p', 'MarkerSize',22, ...
         'MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor','k','LineWidth',1.2);
    text(cfg.Tx1(1)+1.5, cfg.Tx1(2), sprintf('Tx Oeste\n[%g,%g,%.1f]',cfg.Tx1), ...
         'FontSize',11,'FontWeight','bold','VerticalAlignment','middle');
    text(cfg.Tx2(1)-1.5, cfg.Tx2(2), sprintf('Tx Leste\n[%g,%g,%.1f]',cfg.Tx2), ...
         'FontSize',11,'FontWeight','bold','HorizontalAlignment','right','VerticalAlignment','middle');

    % --- RIS de teste (parede Sul) ---
    hR = plot(ris_xy(1), ris_xy(2), 's', 'MarkerSize',16, ...
              'MarkerFaceColor',[0.85 0.30 0.30],'MarkerEdgeColor','k','LineWidth',1.5);
    text(ris_xy(1), ris_xy(2)+6, sprintf('RIS teste [%.0f, %.0f]\nz = %.1f e %.1f m', ...
             ris_xy(1), ris_xy(2), ris_zs(1), ris_zs(2)), ...
         'HorizontalAlignment','center','FontSize',10,'FontWeight','bold','Color',[0.6 0 0], ...
         'BackgroundColor',[1 1 1],'EdgeColor',[0.85 0.30 0.30],'Margin',3);

    % --- rotulos das paredes ---
    text(LX/2, -3, 'Parede Sul (y=0)', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);
    text(LX/2, LY+3, 'Parede Norte (y=50)', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);

    hold off; axis equal;
    xlim([-10 LX+12]); ylim([-8 LY+8]);
    xlabel('x (m)','FontSize',12); ylabel('y (m)','FontSize',12);
    title({'Cenario B - Topologia com RIS de teste (parede Sul)', ...
           sprintf('RIS em x=%.0f, y=%.0f (z testado: %.1f e %.1f m) | obstaculos %d dB | %d usuarios', ...
                   ris_xy(1), ris_xy(2), ris_zs(1), ris_zs(2), cfg.att_obstaculo, size(pontos_rx,1))}, ...
          'FontSize',12);
    legend([hT hU hR], {'Transmissores','Usuarios','RIS de teste'}, ...
           'Location','northeastoutside','FontSize',10);
    grid on; box on;

    fn = 'topologia_cenarioB_ris_teste.png';
    saveas(fig, fn); close(fig);
    fprintf('Gerado: %s\n', fn);
end

% =========================================================================
% topologia_ris_otim (RIS otimas N=64/N=256) -- era plot_topologia_cenarioB_ris_otim.m
% =========================================================================
function fig_topologia_ris_otim()
    thisdir = fileparts(mfilename('fullpath'));
    addpath(thisdir);
    addpath(fileparts(thisdir));   % raiz: SimRIS_v18
    cd(thisdir);

    cfg = cenarioB_config();
    pontos_rx = gera_usuarios_cenarioB(cfg);

    % --- le as melhores pos atuais (N=64 e N=256) da CSV --------------------
    csv_path = fullfile('convergencia_cenarioB', 'progresso_cenarioB.csv');
    T = readmatrix(csv_path);
    Ns = [64, 256];
    cores_ris = {[0.20 0.45 0.85], [0.85 0.30 0.30]};   % N=64 azul, N=256 vermelho
    ris_info = struct();
    for kn = 1:numel(Ns)
        N = Ns(kn);
        sub = T(T(:,1)==N, :);
        if isempty(sub)
            ris_info(kn).N = N; ris_info(kn).has = false; continue;
        end
        [~, ord] = sort(sub(:,2)); sub = sub(ord, :);
        ris_info(kn).N      = N;
        ris_info(kn).has    = true;
        ris_info(kn).it     = sub(end,2);
        ris_info(kn).fitness= sub(end,3);
        ris_info(kn).pos    = sub(end,4:6);
        ris_info(kn).cor    = cores_ris{kn};
        fprintf('N=%d : iter %d | fitness %.3f | best [%.2f, %.2f, %.2f]\n', ...
                N, sub(end,2), sub(end,3), sub(end,4), sub(end,5), sub(end,6));
    end

    LX = cfg.room(1); LY = cfg.room(2);

    fig = figure('Visible','off','Position',[100 100 1100 700]); hold on;
    plot([0 LX LX 0 0], [0 0 LY LY 0], 'k-', 'LineWidth', 2);

    % --- obstaculos ---
    B = cfg.blockers;
    for k = 1:size(B,1)
        x1=B(k,1); x2=B(k,2); y1=B(k,3); y2=B(k,4);
        patch([x1 x2 x2 x1],[y1 y1 y2 y2], [0.95 0.55 0.15], 'FaceAlpha',0.55, ...
              'EdgeColor','k','LineWidth',1.2);
    end

    % --- usuarios ---
    hU = scatter(pontos_rx(:,1), pontos_rx(:,2), 26, [0.20 0.45 0.85], 'filled', ...
                 'MarkerEdgeColor','k','LineWidth',0.3, 'MarkerFaceAlpha',0.55);

    % --- transmissores ---
    hT = plot(cfg.Tx1(1), cfg.Tx1(2), 'p', 'MarkerSize',22, ...
              'MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor','k','LineWidth',1.2);
    plot(cfg.Tx2(1), cfg.Tx2(2), 'p', 'MarkerSize',22, ...
         'MarkerFaceColor',[1 0.6 0],'MarkerEdgeColor','k','LineWidth',1.2);

    % --- marcadores das RIS otimas (sem caixas de anotacao; info na legenda) -
    ris_handles = []; ris_labels = {};
    for kn = 1:numel(ris_info)
        r = ris_info(kn);
        if ~r.has; continue; end
        h = plot(r.pos(1), r.pos(2), 's', 'MarkerSize',16, ...
                 'MarkerFaceColor', r.cor, 'MarkerEdgeColor','k','LineWidth',1.5);
        ris_handles(end+1) = h; %#ok<AGROW>
        ris_labels{end+1} = sprintf('Optimal RIS N=%d  [%.1f, %.1f, %.1f]  fit %.3f', ...
                                    r.N, r.pos(1), r.pos(2), r.pos(3), r.fitness); %#ok<AGROW>
    end

    % --- rotulos das paredes ---
    text(LX/2, -3, 'South Wall (y=0)', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);
    text(LX/2, LY+3, 'North Wall (y=50)', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);

    hold off; axis equal;
    xlim([-10 LX+12]); ylim([-10 LY+8]);
    xlabel('x (m)','FontSize',12); ylabel('y (m)','FontSize',12);
    legend([hT hU ris_handles], [{'Transmitters','Users'} ris_labels], ...
           'Location','northeastoutside','FontSize',10);
    grid on; box on;

    fn = 'topologia_cenarioB_ris_otim.png';
    saveas(fig, fn); close(fig);
    fprintf('Gerado: %s\n', fn);
end

% =========================================================================
% convergencia (FPA N=64 e N=256) -- era plot_convergencia_cenarioB.m
% csv: progresso_cenarioB.csv (pasta convergencia_cenarioB)
% =========================================================================
function fig_convergencia()
    csv = 'progresso_cenarioB.csv';
    T   = readmatrix(csv);                  % N, iter, fitness, x, y, z

    N_lista = [64, 256];
    N_iter  = 300;

    for N = N_lista
        sub = T(T(:,1) == N, :);
        if isempty(sub)
            fprintf('N=%d : no entries in CSV - skipped\n', N);
            continue;
        end
        [~, ord] = sort(sub(:,2));
        sub = sub(ord, :);
        it_real = sub(:,2);
        fb_real = sub(:,3);
        last_it = it_real(end);
        last_fb = fb_real(end);

        fprintf('N=%d : %d iter in CSV (1..%d) | final fbest = %.4f\n', ...
                N, numel(it_real), last_it, last_fb);

        if last_it < N_iter
            it_pad = (last_it+1 : N_iter)';
            fb_pad = repmat(last_fb, numel(it_pad), 1);
            iter_full = [it_real; it_pad];
            fb_full   = [fb_real; fb_pad];
        else
            iter_full = it_real;
            fb_full   = fb_real;
        end

        fig = figure('Visible','off','Position',[100 100 880 540]); hold on;
        plot(iter_full, fb_full, '-o', 'LineWidth', 1.6, 'MarkerSize', 3, ...
             'Color',[0.20 0.45 0.85], 'MarkerFaceColor',[0.20 0.45 0.85]);
        hold off;
        xlabel('Iteration', 'FontSize', 12);
        ylabel('Best Fitness (0.7\cdotmean + 0.3\cdotmin)', 'FontSize', 12);
        xlim([0 N_iter]);
        grid on; box on;

        fn = sprintf('convergencia_cenarioB_N%d.png', N);
        saveas(fig, fn); close(fig);
        fprintf('Generated: %s\n', fn);
    end
end

% =========================================================================
% convergencia_N256 (curva N=256 com padding) -- era convergencia_N256.m
% csv: ../convergencia_cenarioB/progresso_cenarioB.csv
% =========================================================================
function fig_convergencia_N256()
    csv = fullfile('..','convergencia_cenarioB','progresso_cenarioB.csv');
    T   = readmatrix(csv);                  % N, iter, fitness, x, y, z

    N_target  = 256;
    N_iter    = 300;                        % planned target in fpa_cenarioB
    sub = T(T(:,1) == N_target, :);
    [~, ord] = sort(sub(:,2));
    sub = sub(ord, :);
    it_real = sub(:,2);
    fb_real = sub(:,3);
    last_it = it_real(end);
    last_fb = fb_real(end);

    fprintf('N=%d : %d iter in CSV (1..%d) | final fbest = %.4f\n', ...
            N_target, numel(it_real), last_it, last_fb);

    % --- fill iter (last_it+1)..N_iter with last_fb -------------------------
    it_pad = (last_it+1 : N_iter)';
    fb_pad = repmat(last_fb, numel(it_pad), 1);

    iter_full = [it_real; it_pad];
    fb_full   = [fb_real; fb_pad];

    % --- plot ---------------------------------------------------------------
    fig = figure('Visible','off','Position',[100 100 880 540]); hold on;

    plot(iter_full, fb_full, '-o', 'LineWidth', 1.6, 'MarkerSize', 3, ...
         'Color',[0.20 0.45 0.85], 'MarkerFaceColor',[0.20 0.45 0.85]);

    hold off;
    xlabel('Iteration', 'FontSize', 12);
    ylabel('Best Fitness (0.7\cdotmean + 0.3\cdotmin)', 'FontSize', 12);
    xlim([0 N_iter]);
    grid on; box on;

    fn = sprintf('convergencia_cenarioB_N%d.png', N_target);
    saveas(fig, fn); close(fig);
    fprintf('Generated: %s\n', fn);

    % --- export completed CSV -----------------------------------------------
    fn_csv = sprintf('convergencia_cenarioB_N%d_completo.csv', N_target);
    fid = fopen(fn_csv, 'w');
    fprintf(fid, 'iter,fbest,source\n');
    for k = 1:numel(iter_full)
        if iter_full(k) <= last_it
            src = 'real';
        else
            src = 'padded';
        end
        fprintf(fid, '%d,%.6f,%s\n', iter_full(k), fb_full(k), src);
    end
    fclose(fid);
    fprintf('Generated: %s (%d rows)\n', fn_csv, numel(iter_full));
end

% =========================================================================
% boxplot_otim_parcial (parcial da otimizacao) -- era boxplot_otim_parcial.m
% =========================================================================
function fig_boxplot_otim_parcial()
    rng(1);
    thisdir = fileparts(mfilename('fullpath'));
    addpath(thisdir);
    addpath(fileparts(thisdir));   % root: SimRIS_v18
    cd(thisdir);

    cfg = cenarioB_config();
    pos_Rx = gera_usuarios_cenarioB(cfg);

    % --- reused baseline (from teste_ris_alta_cenarioB) ---------------------
    S = load('teste_ris_alta_cenarioB.mat', 'grupos', 'rotulos');
    taxas_base = S.grupos{1};       % 96x1
    fprintf('Baseline (reused): median %.3f | mean %.3f | min %.3f\n', ...
            median(taxas_base), mean(taxas_base), min(taxas_base));

    % --- read progress CSV and pick the best pos per N ----------------------
    csv_path = fullfile('convergencia_cenarioB', 'progresso_cenarioB.csv');
    T = readmatrix(csv_path);       % N, iter, fitness, x, y, z
    Ns_target = [64, 256];

    grupos = {taxas_base};
    rotulos = {'Baseline (no RIS)'};
    posicoes = {'No RIS'};
    its = {[]};

    for N = Ns_target
        sub = T(T(:,1) == N, :);
        if isempty(sub)
            fprintf('N=%d : no iterations yet - skipped from boxplot\n', N);
            continue;
        end
        [~, ord] = sort(sub(:,2));
        sub = sub(ord, :);
        it_max = sub(end, 2);
        best = sub(end, 4:6);
        fbest = sub(end, 3);
        sc = 1;                     % South/North side walls -> Scenario 1
        fprintf('N=%d : iter %d | best [%.2f, %.2f, %.2f] | fitness %.3f\n', ...
                N, it_max, best(1), best(2), best(3), fbest);
        [~, taxas] = calculamediavazao_cenarioB(best, sc, N, pos_Rx, cfg);
        grupos{end+1}   = taxas(:); %#ok<AGROW>
        rotulos{end+1}  = sprintf('RIS N=%d (it %d)', N, it_max); %#ok<AGROW>
        posicoes{end+1} = sprintf('[%.1f, %.1f, %.1f]', best(1),best(2),best(3)); %#ok<AGROW>
        its{end+1}      = it_max; %#ok<AGROW>
        fprintf('  rates: median %.3f | mean %.3f | min %.3f\n', ...
                median(taxas), mean(taxas), min(taxas));
    end

    % --- Manual boxplot (no Statistics Toolbox) -----------------------------
    fig = figure('Visible','off','Position',[100 100 820 580]); hold on;
    cores = [0.70 0.70 0.75; 0.20 0.45 0.85; 0.85 0.30 0.30];
    larg = 0.5; ymax = 0; topo = zeros(1, numel(grupos));
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
        if isempty(outl); topo(g) = hi; else; topo(g) = max(hi, max(outl)); end
    end
    % position annotation above each box
    for g = 1:numel(grupos)
        text(g, topo(g) + ymax*0.045, posicoes{g}, 'HorizontalAlignment','center', ...
             'FontSize', 9, 'FontWeight','bold', 'Color',[0.15 0.15 0.5]);
    end
    hold off;
    set(gca,'XTick',1:numel(rotulos),'XTickLabel',rotulos,'FontSize',11);
    xlim([0.4 numel(rotulos)+0.6]); ylim([0 max(ymax*1.18, 0.1)]);
    ylabel('Spectral Efficiency (bps/Hz)','FontSize',12);
    legend(hM, {'Mean'}, 'Location','northeast','FontSize',10);
    grid on; box on;

    fn = 'boxplot_cenarioB_otim_parcial.png';
    saveas(fig, fn); close(fig);
    fprintf('Generated: %s\n', fn);
end

% =========================================================================
% boxplot_sinr (SE multi-Tx SINR, recomputa) -- era boxplot_sinr_cenarioB.m
% =========================================================================
function fig_boxplot_sinr()
    rng(1);
    thisdir = fileparts(mfilename('fullpath'));
    addpath(thisdir);
    addpath(fileparts(thisdir));   % root: SimRIS_v18
    cd(thisdir);

    cfg = cenarioB_config();
    pos_Rx = gera_usuarios_cenarioB(cfg);

    % --- read progress CSV and pick the best pos per N ----------------------
    csv_path = fullfile('convergencia_cenarioB', 'progresso_cenarioB.csv');
    T = readmatrix(csv_path);       % N, iter, fitness, x, y, z
    Ns_target = [64, 256];
    sc = 1;                         % South/North side walls -> Scenario 1

    best_by_N = containers.Map('KeyType','double','ValueType','any');
    for N = Ns_target
        sub = T(T(:,1) == N, :);
        if isempty(sub)
            error('N=%d has no entries in CSV.', N);
        end
        % best fitness (on tie, last iteration)
        [fbest, ~] = max(sub(:,3));
        cands = sub(sub(:,3) == fbest, :);
        [~, k] = max(cands(:,2));
        chosen = cands(k, :);
        best_by_N(N) = struct('it', chosen(2), 'pos', chosen(4:6), 'fit', fbest);
        fprintf('N=%3d : iter %3d | best [%.2f, %.2f, %.2f] | fitness %.4f\n', ...
                N, chosen(2), chosen(4), chosen(5), chosen(6), fbest);
    end

    % --- recompute baseline + RIS for each N --------------------------------
    grupos    = {};
    rotulos   = {};
    posicoes  = {};
    taxas_base = [];

    for N = Ns_target
        info = best_by_N(N);
        fprintf('\n>>> computing N=%d (pos [%.2f, %.2f, %.2f]) ...\n', ...
                N, info.pos(1), info.pos(2), info.pos(3));
        t0 = tic;
        [~, taxas, taxas_bl] = calculamediavazao_cenarioB(info.pos, sc, N, pos_Rx, cfg);
        fprintf('    done in %.1fs | RIS: median %.3f mean %.3f min %.3f\n', ...
                toc(t0), median(taxas), mean(taxas), min(taxas));
        fprintf('    BL: median %.3f mean %.3f min %.3f\n', ...
                median(taxas_bl), mean(taxas_bl), min(taxas_bl));
        if isempty(taxas_base)
            taxas_base = taxas_bl(:);
        end
        grupos{end+1}   = taxas(:); %#ok<AGROW>
        rotulos{end+1}  = sprintf('RIS N=%d', N); %#ok<AGROW>
        posicoes{end+1} = sprintf('[%.1f, %.1f, %.1f]', info.pos(1), info.pos(2), info.pos(3)); %#ok<AGROW>
    end

    % --- prepend baseline ----------------------------------------------------
    grupos    = [{taxas_base},   grupos];
    rotulos   = [{'Baseline (no RIS)'}, rotulos];
    posicoes  = [{'No RIS'},    posicoes];

    % --- Manual boxplot (no Statistics Toolbox) -----------------------------
    fig = figure('Visible','off','Position',[100 100 900 600]); hold on;
    cores = [0.70 0.70 0.75; 0.20 0.45 0.85; 0.85 0.30 0.30];
    larg = 0.5; ymax = 0; topo = zeros(1, numel(grupos));
    hM = [];
    for g = 1:numel(grupos)
        v = sort(grupos{g});
        q1 = pct(v,25); med = pct(v,50); q3 = pct(v,75); iqr = q3-q1;
        lo = max(min(v), q1-1.5*iqr); hi = min(max(v), q3+1.5*iqr);
        outl = v(v<lo | v>hi);
        cor = cores(mod(g-1,size(cores,1))+1,:);
        xL = g-larg/2; xR = g+larg/2;
        fill([xL xR xR xL],[q1 q1 q3 q3], cor, 'FaceAlpha',0.6,'EdgeColor','k','LineWidth',1);
        plot([xL xR],[med med],'k-','LineWidth',2);
        hMg = plot(g, mean(grupos{g}),'d','MarkerFaceColor',[0 0.5 0],'MarkerEdgeColor','k','MarkerSize',7);
        if g == 1; hM = hMg; end
        plot([g g],[q3 hi],'k--','LineWidth',1); plot([g g],[lo q1],'k--','LineWidth',1);
        plot([g-0.15 g+0.15],[hi hi],'k-'); plot([g-0.15 g+0.15],[lo lo],'k-');
        if ~isempty(outl); plot(g*ones(size(outl)), outl,'r+','MarkerSize',5); end
        ymax = max(ymax, max(v));
        if isempty(outl); topo(g) = hi; else; topo(g) = max(hi, max(outl)); end
    end
    hold off;
    set(gca,'XTick',1:numel(rotulos),'XTickLabel',rotulos,'FontSize',11);
    xlim([0.4 numel(rotulos)+0.6]); ylim([0 max(ymax*1.18, 0.1)]);
    ylabel('SINR','FontSize',12);
    legend(hM, {'Mean'}, 'Location','northeast','FontSize',10);
    grid on; box on;

    fn_png = 'boxplot_sinr_cenarioB.png';
    fn_mat = 'boxplot_sinr_cenarioB.mat';
    saveas(fig, fn_png); close(fig);
    save(fn_mat, 'grupos', 'rotulos', 'posicoes', 'best_by_N', 'cfg');
    fprintf('\nGenerated: %s | %s\n', fn_png, fn_mat);
end

% =========================================================================
% boxplot_mc (SE from Monte Carlo CSV) -- era boxplot_mc_cenarioB.m
% =========================================================================
function fig_boxplot_mc()
    thisdir = fileparts(mfilename('fullpath'));
    cd(thisdir);

    csv = 'mc_se_ues_cenarioB.csv';
    T = readtable(csv);
    fprintf('Read %s : %d rows\n', csv, height(T));

    grupos   = {T.baseline, T.RIS_N64, T.RIS_N256};
    rotulos  = {'Baseline (no RIS)', 'RIS N=64', 'RIS N=256'};
    posicoes = {'No RIS', '[59.2, 50.0, 3.0]', '[56.2, 50.0, 3.0]'};
    nG = numel(grupos);

    % --- summary ------------------------------------------------------------
    fprintf('\nSE summary (bps/Hz):\n');
    fprintf('%-22s %8s %8s %8s %8s %8s %8s\n', ...
            'Group','min','P5','P25','median','mean','max');
    for g = 1:nG
        v = sort(grupos{g}(~isnan(grupos{g})));
        fprintf('%-22s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', rotulos{g}, ...
            v(1), pct(v,5), pct(v,25), pct(v,50), mean(v), v(end));
    end

    % --- manual boxplot full -------------------------------------------------
    plot_box(grupos, rotulos, posicoes, 'boxplot_mc_cenarioB.png', []);
    fprintf('Generated: boxplot_mc_cenarioB.png\n');

    % --- manual boxplot zoom 0..4 -------------------------------------------
    plot_box(grupos, rotulos, posicoes, 'boxplot_mc_cenarioB_zoom.png', [0 4]);
    fprintf('Generated: boxplot_mc_cenarioB_zoom.png\n');
end

% =========================================================================
% cdf_mc (ECDF from Monte Carlo CSV) -- era cdf_mc_cenarioB.m
% =========================================================================
function fig_cdf_mc()
    thisdir = fileparts(mfilename('fullpath'));
    cd(thisdir);

    csv = 'mc_se_ues_cenarioB.csv';
    T = readtable(csv);             % cols: realization, UE, baseline, RIS_N64, RIS_N256
    fprintf('Read %s : %d rows\n', csv, height(T));

    grupos  = {T.baseline, T.RIS_N64, T.RIS_N256};
    rotulos = {'Baseline (no RIS)', 'RIS N=64', 'RIS N=256'};
    nG = numel(grupos);

    % --- summary ------------------------------------------------------------
    fprintf('\nSE summary (bps/Hz):\n');
    fprintf('%-22s %8s %8s %8s %8s %8s %8s\n', ...
            'Group','min','P5','P25','median','mean','max');
    for g = 1:nG
        v = sort(grupos{g}(~isnan(grupos{g})));
        fprintf('%-22s %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n', rotulos{g}, ...
            v(1), pct(v,5), pct(v,25), pct(v,50), mean(v), v(end));
    end

    cores = [0.30 0.30 0.35;
             0.20 0.45 0.85;
             0.85 0.30 0.30];
    estilos  = {'-.','--','-'};
    larguras = [1.8, 2.0, 2.4];

    % --- ECDF full -----------------------------------------------------------
    fig = figure('Visible','off','Position',[100 100 950 620]); hold on;
    h = gobjects(1,nG);
    for g = 1:nG
        v = sort(grupos{g}(~isnan(grupos{g})));
        n = numel(v);
        x = [v(1); v];
        F = (0:n)' / n;
        cor = cores(g,:);
        h(g) = stairs(x, F, estilos{g}, 'Color', cor, 'LineWidth', larguras(g));
    end
    plot(xlim, [0.5 0.5],  ':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    plot(xlim, [0.05 0.05],':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    hold off;
    xlabel('SE (bps/Hz) - multi-Tx SINR model','FontSize',12);
    ylabel('Empirical CDF  F(x) = P(SE \leq x)','FontSize',12);
    legend(h, rotulos, 'Location','southeast', 'FontSize',11);
    ylim([0 1]); grid on; box on;
    saveas(fig, 'cdf_mc_cenarioB.png'); close(fig);
    fprintf('Generated: cdf_mc_cenarioB.png\n');

    % --- ECDF zoom 0..3 -----------------------------------------------------
    fig2 = figure('Visible','off','Position',[100 100 950 620]); hold on;
    h2 = gobjects(1,nG);
    for g = 1:nG
        v = sort(grupos{g}(~isnan(grupos{g})));
        n = numel(v);
        x = [v(1); v];
        F = (0:n)' / n;
        cor = cores(g,:);
        h2(g) = stairs(x, F, estilos{g}, 'Color', cor, 'LineWidth', larguras(g));
    end
    xlim([0 3]); ylim([0 1]);
    plot(xlim, [0.5 0.5],  ':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    plot(xlim, [0.05 0.05],':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    hold off;
    xlabel('SE (bps/Hz) - multi-Tx SINR model','FontSize',12);
    ylabel('Empirical CDF  F(x) = P(SE \leq x)','FontSize',12);
    legend(h2, rotulos, 'Location','southeast', 'FontSize',11);
    grid on; box on;
    saveas(fig2, 'cdf_mc_cenarioB_zoom.png'); close(fig2);
    fprintf('Generated: cdf_mc_cenarioB_zoom.png\n');

    % --- ECDF log-X to visualize the low tail -------------------------------
    fig3 = figure('Visible','off','Position',[100 100 950 620]); hold on;
    h3 = gobjects(1,nG);
    for g = 1:nG
        v = sort(grupos{g}(~isnan(grupos{g})));
        v = v(v > 0);                    % log does not accept 0
        n = numel(v);
        x = [v(1); v];
        F = (0:n)' / n;
        cor = cores(g,:);
        h3(g) = stairs(x, F, estilos{g}, 'Color', cor, 'LineWidth', larguras(g));
    end
    set(gca,'XScale','log');
    plot(xlim, [0.5 0.5], ':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    plot(xlim, [0.05 0.05],':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    hold off;
    xlabel('SE (bps/Hz) - log scale','FontSize',12);
    ylabel('Empirical CDF  F(x) = P(SE \leq x)','FontSize',12);
    legend(h3, rotulos, 'Location','southeast', 'FontSize',11);
    ylim([0 1]); grid on; box on;
    saveas(fig3, 'cdf_mc_cenarioB_logx.png'); close(fig3);
    fprintf('Generated: cdf_mc_cenarioB_logx.png\n');
end

% =========================================================================
% cdf_sinr (ECDF per-UE from boxplot_sinr mat) -- era cdf_sinr_cenarioB.m
% =========================================================================
function fig_cdf_sinr()
    thisdir = fileparts(mfilename('fullpath'));
    cd(thisdir);

    S = load('boxplot_sinr_cenarioB.mat', 'grupos', 'rotulos', 'posicoes');
    grupos   = S.grupos;
    rotulos  = S.rotulos;
    nG = numel(grupos);
    nU = numel(grupos{1});
    fprintf('UEs per group: %d\n', nU);

    % --- CSV table: SE of each UE ------------------------------------------
    H = {'UE'};
    for g = 1:nG
        H{end+1} = matlab.lang.makeValidName(rotulos{g}); %#ok<AGROW>
    end
    M = (1:nU)';
    for g = 1:nG
        M = [M, grupos{g}(:)]; %#ok<AGROW>
    end
    fid = fopen('se_por_ue.csv','w');
    fprintf(fid, '%s', H{1}); for k=2:numel(H); fprintf(fid, ',%s', H{k}); end
    fprintf(fid, '\n');
    for i=1:nU
        fprintf(fid, '%d', M(i,1));
        for k=2:size(M,2); fprintf(fid, ',%.6f', M(i,k)); end
        fprintf(fid, '\n');
    end
    fclose(fid);
    fprintf('Generated: se_por_ue.csv (%d UEs x %d groups)\n', nU, nG);

    % --- statistical summary ------------------------------------------------
    fprintf('\nSE summary (bps/Hz):\n');
    fprintf('%-25s %8s %8s %8s %8s %8s %8s\n', ...
            'Group','min','P5','P25','median','mean','max');
    for g = 1:nG
        v = sort(grupos{g}(:));
        fprintf('%-25s %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', rotulos{g}, ...
            v(1), pct(v,5), pct(v,25), pct(v,50), mean(v), v(end));
    end

    % --- ECDF ---------------------------------------------------------------
    cores = [0.30 0.30 0.35;    % baseline
             0.20 0.45 0.85;    % N=64
             0.85 0.30 0.30];   % N=256
    estilos  = {'-.','--','-'};
    larguras = [1.8, 2.0, 2.4];
    markers  = {'o','s','^'};

    fig = figure('Visible','off','Position',[100 100 950 620]); hold on;
    h = gobjects(1,nG);
    for g = 1:nG
        v = sort(grupos{g}(:));
        n = numel(v);
        x = [v(1); v];          % anchor at min for y=0
        F = (0:n)' / n;         % 0, 1/n, 2/n, ..., 1
        cor = cores(mod(g-1,size(cores,1))+1,:);
        h(g) = stairs(x, F, estilos{mod(g-1,numel(estilos))+1}, ...
            'Color', cor, 'LineWidth', larguras(mod(g-1,numel(larguras))+1));
        % markers at each UE (1 point per UE)
        plot(v, F(2:end), markers{mod(g-1,numel(markers))+1}, ...
            'MarkerEdgeColor', cor, 'MarkerFaceColor','none', 'MarkerSize',4);
    end
    % horizontal guides
    plot(xlim, [0.5 0.5], ':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    plot(xlim, [0.05 0.05], ':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    hold off;
    xlabel('SE per user (bps/Hz) - SINR model','FontSize',12);
    ylabel('Empirical CDF  F(x) = P(SE \leq x)','FontSize',12);
    legend(h, rotulos, 'Location','southeast', 'FontSize',11);
    ylim([0 1]); grid on; box on;
    fn = 'cdf_sinr_cenarioB.png';
    saveas(fig, fn); close(fig);
    fprintf('\nGenerated: %s\n', fn);

    % --- zoom 0..2.5 --------------------------------------------------------
    fig2 = figure('Visible','off','Position',[100 100 950 620]); hold on;
    h2 = gobjects(1,nG);
    for g = 1:nG
        v = sort(grupos{g}(:));
        n = numel(v);
        x = [v(1); v];
        F = (0:n)' / n;
        cor = cores(mod(g-1,size(cores,1))+1,:);
        h2(g) = stairs(x, F, estilos{mod(g-1,numel(estilos))+1}, ...
            'Color', cor, 'LineWidth', larguras(mod(g-1,numel(larguras))+1));
        plot(v, F(2:end), markers{mod(g-1,numel(markers))+1}, ...
            'MarkerEdgeColor', cor, 'MarkerFaceColor','none', 'MarkerSize',5);
    end
    xlim([0 2.5]); ylim([0 1]);
    plot(xlim, [0.5 0.5], ':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    plot(xlim, [0.05 0.05], ':', 'Color',[0.5 0.5 0.5], 'LineWidth',0.8);
    hold off;
    xlabel('SE per user (bps/Hz) - SINR model','FontSize',12);
    ylabel('Empirical CDF  F(x) = P(SE \leq x)','FontSize',12);
    legend(h2, rotulos, 'Location','southeast', 'FontSize',11);
    grid on; box on;
    fn2 = 'cdf_sinr_cenarioB_zoom.png';
    saveas(fig2, fn2); close(fig2);
    fprintf('Generated: %s\n', fn2);
end

% =========================================================================
% HELPERS COMPARTILHADOS (antes duplicados em varios scripts)
% =========================================================================

% Boxplot manual (sem Statistics Toolbox). Usado por boxplot_mc.
function plot_box(grupos, rotulos, posicoes, fn, ylim_custom)
    cores = [0.70 0.70 0.75;
             0.20 0.45 0.85;
             0.85 0.30 0.30];
    fig = figure('Visible','off','Position',[100 100 900 620]); hold on;
    larg = 0.5; topo = zeros(1, numel(grupos)); ymax = 0;
    hM = [];
    for g = 1:numel(grupos)
        v = sort(grupos{g}(~isnan(grupos{g})));
        q1 = pct(v,25); med = pct(v,50); q3 = pct(v,75); iqr = q3-q1;
        lo = max(min(v), q1-1.5*iqr); hi = min(max(v), q3+1.5*iqr);
        outl = v(v<lo | v>hi);
        cor = cores(mod(g-1,size(cores,1))+1,:);
        xL = g-larg/2; xR = g+larg/2;
        fill([xL xR xR xL],[q1 q1 q3 q3], cor, ...
             'FaceAlpha',0.6,'EdgeColor','k','LineWidth',1);
        plot([xL xR],[med med],'k-','LineWidth',2);
        hMg = plot(g, mean(v),'d','MarkerFaceColor',[0 0.5 0], ...
                   'MarkerEdgeColor','k','MarkerSize',7);
        if g == 1; hM = hMg; end
        plot([g g],[q3 hi],'k--','LineWidth',1);
        plot([g g],[lo q1],'k--','LineWidth',1);
        plot([g-0.15 g+0.15],[hi hi],'k-');
        plot([g-0.15 g+0.15],[lo lo],'k-');
        if ~isempty(outl)
            plot(g*ones(size(outl)), outl,'r+','MarkerSize',4);
        end
        ymax = max(ymax, max(v));
        if isempty(outl); topo(g) = hi; else; topo(g) = max(hi, max(outl)); end
    end
    if isempty(ylim_custom)
        yl = [0 max(ymax*1.18, 0.1)];
    else
        yl = ylim_custom;
    end
    ylim(yl);
    % annotations above boxes (only if they fit in ylim)
    for g = 1:numel(grupos)
        ypos = min(topo(g) + (yl(2)-yl(1))*0.045, yl(2)*0.97);
        text(g, ypos, posicoes{g}, 'HorizontalAlignment','center', ...
             'FontSize',9,'FontWeight','bold','Color',[0.15 0.15 0.5]);
    end
    hold off;
    set(gca,'XTick',1:numel(rotulos),'XTickLabel',rotulos,'FontSize',11);
    xlim([0.4 numel(rotulos)+0.6]);
    ylabel('SE per sample (bps/Hz) - SINR model','FontSize',12);
    legend(hM, {'Mean'}, 'Location','northeast','FontSize',10);
    grid on; box on;
    saveas(fig, fn); close(fig);
end

% Percentil (interpolacao linear), sem Statistics Toolbox.
function q = pct(v, p)
    n = numel(v);
    if n==1; q=v(1); return; end
    pos = (p/100)*n + 0.5; pos = min(max(pos,1),n);
    lo = floor(pos); hi = ceil(pos);
    if lo==hi; q=v(lo); else; q = v(lo)+(pos-lo)*(v(hi)-v(lo)); end
end
