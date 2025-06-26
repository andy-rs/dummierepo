%% Métodos de selección

% Torneo
function seleccionados = seleccionPorTorneo(poblacion, fitness, k)
    len = size(poblacion, 1);
    seleccionados = zeros(size(poblacion));
    
    for i = 1:len
        indices = randsample(len, k);
        [~, best_idx] = max(fitness(indices));
        seleccionados(i,:) = poblacion(indices(best_idx), :);
    end
end

% Ruleta

function seleccionados = seleccionPorRuleta(poblacion, fitness)
    len = size(poblacion, 1);
    probs = fitness / sum(fitness);
    indices_seleccionados = randsample(1:len, len, true, probs);
    seleccionados = poblacion(indices_seleccionados, :);
end

% Ranking
function seleccionados = seleccionPorRanking(poblacion, fitness)
    len = size(poblacion, 1);
    [~, orden] = sort(fitness, 'ascend');
    rangos = 1:len;
    probs = rangos / sum(rangos);
    indices_seleccionados = randsample(orden, len, true, probs);
    seleccionados = poblacion(indices_seleccionados, :); 
end

% Elite
function seleccionados = seleccionPorElitismo(poblacion, fitness, num_seleccionados)
    len = size(poblacion, 1);
    [~, orden] = sort(fitness, 'descend');
    indices_seleccionados = orden(1:num_seleccionados);
    seleccionados = poblacion(indices_seleccionados, :);
end
