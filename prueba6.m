%% Códigos de práctica para algoritmos genéticos

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

%% Métodos de cruce

% Cruce de un punto

function hijos = cruceDeUnPunto(poblacion)
    hijos = zeros(size(poblacion));
    len = size(poblacion, 1);

    num_genes = size(poblacion, 2);

    for i = 1:2:len
        p1 = poblacion(i, :);
        p2 = poblacion(i+1, :);
        punto = randi([1 num_genes-1]);
        hijos(i, :) = [p1(1:punto) p2(punto+1:end)];
        hijos(i+1, :) = [p2(1:punto) p1(punto+1:end)];
    end
end

% Cruce de dos puntos

function hijos = cruceDeDosPuntos(poblacion)
    len = size(poblacion, 1);
    num_genes = size(poblacion,2);
    
    hijos = zeros(len, num_genes);

    for i = 1:2:len
        p1 = poblacion(i, :);
        p2 = poblacion(i+1, :);

        puntos = sort(randsample(1:num_genes-1, 2));
        
        hijos(i, :) = [p1(1:puntos(1)) p2(puntos(1)+1:puntos(2)) p1(puntos(2)+1:end)];
        hijos(i+1, :) = [p2(1:puntos(1)) p1(puntos(1)+1:puntos(2)) p2(puntos(2)+1:end)];
    end
end

% Cruce uniforme

function hijos = cruceUniforme(poblacion, probabilidad)
    [len, num_genes] = size(poblacion);
    
    hijos = zeros(len, num_genes);

    for i = 1:2:len
        p1 = poblacion(i, :);
        p2 = poblacion(i+1, :);

        mascara = rand(1, num_genes) < probabilidad;
        
        hijos(i, :) = p1;
        hijos(i+1, :) = p2;

        hijos(i, mascara) = p2(mascara);
        hijos(i+1, mascara) = p1(mascara);
    end
end

% Cruce aritmético

function hijos = cruceAritmetico(seleccionados)
    len = size(seleccionados, 1);
    num_genes = size(seleccionados, 2);

    hijos = zeros(len, num_genes);
    
    for i = 1:2:len
        p1 = seleccionados(i, :);
        p2 = seleccionados(i+1, :);

        alfa = rand();

        hijos(i, :) = alfa * p1 + (1-alfa) * p2;
        hijos(i+1, :) = (1-alfa) * p1 + alfa * p2;
    end
end

%% Métodos de mutación

% Mutación flipBit
function mutados = mutacionFlipBit(hijos, probabilidad)
    [len, num_genes] = size(hijos);
    mascara = rand(len, num_genes) < probabilidad;

    mutados = hijos;
    mutados(mascara) = 1 - mutados(mascara);
end

% Mutacion Uniforme
function mutados = mutacionUniforme(hijos, probabilidad, rango_min, rango_max)
    [len, num_genes] = size(hijos);

    mascara = rand(len, num_genes) < probabilidad;

    mutados = hijos;
    alfa = rand(len, num_genes);
    mutados(mascara) = rango_min + (rango_max - rango_min) * alfa(mascara);
end

% Mutación por intercambio
function mutados = mutacionPorIntercambio(hijosm prob)
    [len, num_genes] = size(hijos);

    mutados = hijos; 
    for i = 1:len
        if rand() < prob
            puntos = randsample(1:num_genes, 2);
            temp = mutados(i, puntos(1));
            mutados(i, puntos(1)) = mutados(i, puntos(2));
            mutados(i, puntos(2)) = temp; 
        end
    end
end

% Mutación Gaussiana
function mutados = mutacionGaussiana(seleccionados, sigma, prob)
    [len, num_genes] = size(seleccionados);
    mask = rand(len, num_genes) < prob;

    mutados = seleccionados;
    ruido = sigma * randn(len, num_genes);
    mutados(mask) = mutados(mask) + ruido(mask);
end