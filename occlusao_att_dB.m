function att_dB = occlusao_att_dB(P1, P2, blockers)
% Atenuacao extra (dB) sofrida pelo segmento P1->P2 ao atravessar caixas
% (obstaculos / paredes grossas). Cada bloqueador e uma AABB 3D:
%   blockers(k,:) = [xmin xmax ymin ymax zmin zmax att_dB]
% Teste de intersecao segmento-AABB pelo metodo dos slabs. Se o raio
% atravessa a caixa dentro do trecho [0,1], soma-se att_dB do bloqueador.
%
% P1, P2 : pontos 3D (linha 1x3), em coordenadas REAIS da sala.
% att_dB : perda total de penetracao acumulada (>= 0).

    att_dB = 0;
    if isempty(blockers); return; end
    d = P2 - P1;
    EPS = 1e-9;

    for k = 1:size(blockers,1)
        bmin = blockers(k, [1 3 5]);
        bmax = blockers(k, [2 4 6]);
        tmin = 0; tmax = 1; hit = true;
        for ax = 1:3
            if abs(d(ax)) < EPS
                if P1(ax) < bmin(ax) || P1(ax) > bmax(ax)
                    hit = false; break;   % paralelo e fora do slab
                end
            else
                t1 = (bmin(ax) - P1(ax)) / d(ax);
                t2 = (bmax(ax) - P1(ax)) / d(ax);
                if t1 > t2; tmp = t1; t1 = t2; t2 = tmp; end
                tmin = max(tmin, t1);
                tmax = min(tmax, t2);
                if tmin > tmax; hit = false; break; end
            end
        end
        % conta so travessia real (interseccao de comprimento positivo);
        % toque de quina/aresta (tmin==tmax) NAO bloqueia.
        if hit && (tmax - tmin) > 1e-9
            att_dB = att_dB + blockers(k, 7);
        end
    end
end
