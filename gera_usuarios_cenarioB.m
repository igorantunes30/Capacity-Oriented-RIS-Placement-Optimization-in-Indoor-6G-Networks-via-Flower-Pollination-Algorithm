function pos_Rx = gera_usuarios_cenarioB(cfg)
% Sorteia cfg.num_users usuarios uniformemente na sala (semente fixa),
% rejeitando pontos dentro da pegada (x-y) dos bloqueadores, perto demais
% das paredes (cfg.margem) ou dos transmissores (cfg.dist_min_tx).

    rng(cfg.seed_users);
    Lx = cfg.room(1); Ly = cfg.room(2);
    m  = cfg.margem;
    z  = cfg.user_z;

    pos_Rx = zeros(cfg.num_users, 3);
    n = 0; tentativas = 0; max_tent = cfg.num_users * 500;

    while n < cfg.num_users && tentativas < max_tent
        tentativas = tentativas + 1;
        x = m + rand*(Lx - 2*m);
        y = m + rand*(Ly - 2*m);

        % rejeita se cair dentro da pegada x-y de algum bloqueador
        dentro = false;
        for k = 1:size(cfg.blockers,1)
            b = cfg.blockers(k,:);
            if x >= b(1) && x <= b(2) && y >= b(3) && y <= b(4)
                dentro = true; break;
            end
        end
        if dentro; continue; end

        % rejeita se perto demais de algum Tx (no plano)
        if norm([x y] - cfg.Tx1(1:2)) < cfg.dist_min_tx; continue; end
        if norm([x y] - cfg.Tx2(1:2)) < cfg.dist_min_tx; continue; end

        n = n + 1;
        pos_Rx(n,:) = [x, y, z];
    end

    if n < cfg.num_users
        error('gera_usuarios_cenarioB: so foram colocados %d/%d usuarios.', ...
              n, cfg.num_users);
    end
end
