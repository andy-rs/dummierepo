% Solución (Mochila)

% Parámetros del problema
pesos = [10, 20, 30, 15, 25, 5, 35, 12, 22, 18];
valores = [60, 100, 120, 70, 90, 30, 150, 50, 80, 110];
capacidad_maxima = 100;
num_items = length(pesos);

% Parámetros del algoritmo genético
num_generaciones = 200;
tam_poblacion = 100;
prob_cruce = 0.8;
prob_mutacion = 0.02;

maximos = zeros(num_generaciones, 1);
mejores_individuos = zeros(num_generaciones, num_items);
% Inicialización de la población
poblacion = randi([0 1], tam_poblacion, num_items);

for gen = 1:num_generaciones
    % evaluación de aptitud
    aptitudes = evaluarAptitud(poblacion, pesos, valores);
    % mejor aptitud
    maximos(gen)= max(aptitudes);
    [~, idx] = max(aptitudes);
    mejores_individuos(gen, :) = poblacion(idx, :);
    % selección
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 6);
    % cruce 
    hijos = cruceDeUnPunto(poblacion, prob_cruce);
    % mutacion
    mutados = mutacionFlipBit(hijos, prob_mutacion);
    % Reemplazar poblacion original
    poblacion = mutados; 
end

function aptitudes = evaluarAptitud(poblacion, pesos, valores)
    [len, num_genes] = size(poblacion);
    aptitudes = zeros(len, 1);
    
    for i = 1:len
        mask = poblacion(i,:) == 1;
        pesos_individuo = pesos(mask);
        valor_individuo = sum(valores(mask));
        peso_total = sum(pesos_individuo);
        if peso_total > 100
            aptitudes(i) = 0;
        else
            aptitudes(i) = valor_individuo;
        end
    end
end

function seleccionados = seleccionPorTorneo(poblacion, fitness, k)
    len = size(poblacion, 1);
    seleccionados = zeros(size(poblacion));
    
    for i = 1:len
        indices = randsample(len, k);
        [~, best_idx] = max(fitness(indices));
        seleccionados(i,:) = poblacion(indices(best_idx), :);
    end
end

function hijos = cruceDeUnPunto(poblacion, prob)
    hijos = zeros(size(poblacion));
    len = size(poblacion, 1);

    num_genes = size(poblacion, 2);

    for i = 1:2:len
        p1 = poblacion(i, :);
        p2 = poblacion(i+1, :);
        if rand < prob
            punto = randi([1 num_genes-1]);
            hijos(i, :) = [p1(1:punto) p2(punto+1:end)];
            hijos(i+1, :) = [p2(1:punto) p1(punto+1:end)];
        else
            hijos(i, :) = p1;
            hijos(i+1, :) = p2;
        end
    end
end

function mutados = mutacionFlipBit(hijos, probabilidad)
    [len, num_genes] = size(hijos);
    mascara = rand(len, num_genes) < probabilidad;

    mutados = hijos;
    mutados(mascara) = 1 - mutados(mascara);
end

% Gráficas
figure(1)
plot(1:num_generaciones, maximos)
xlabel('Generación')
ylabel('Mejor aptitud')
title('Mejor aptitud de acuerdo a la generacion')

[mejor_aptitud, idx] = max(evaluarAptitud(poblacion, pesos, valores));
mejor_individuo = poblacion(idx, :);
fprintf('Mejor aptitud en la generación que terminó el algoritmo: %d\n', mejor_aptitud);
fprintf('Items seleccionados: %s\n\n\n', mat2str(find(mejor_individuo)));

mejor_aptitud_global = max(maximos);
[~, idy] = max(maximos);
mejor_individuo_global = mejores_individuos(idy, :);
fprintf('La mejor aptitud obtenida es de %d\n', mejor_aptitud_global);
fprintf('Con los items seleccionados: %s\n', mat2str(find(mejor_individuo_global)))

%% Parámetros del problema Himmelblau
rango_min = -5;
rango_max = 5;
n_dim = 2; % x, y

% Parámetros del algoritmo genético
tam_poblacion = 100;
num_generaciones = 200;
prob_cruce = 0.8;
prob_mutacion = 0.2;

% Inicialización
poblacion = rango_min + (rango_max - rango_min) * rand(tam_poblacion, n_dim);
maximos = zeros(num_generaciones, 1);
mejores_individuos = zeros(num_generaciones, n_dim);

for gen = 1:num_generaciones
    % Evaluar aptitud
    aptitudes = evaluarAptitudHimmelblau(poblacion);
    costos = (1 ./ aptitudes) - 1;

    % Guardar mejor
    [~, idx] = max(aptitudes);
    maximos(gen) = costos(idx);
    mejores_individuos(gen, :) = poblacion(idx, :);
    mejor_individuo = poblacion(idx, :);

    % Selección
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 10);

    % Cruce aritmético
    hijos = cruceAritmetico(seleccionados, prob_cruce);

    % Mutación gaussiana
    mutados = mutacionGaussiana(hijos, prob_mutacion, 0.1);

    % Elitismo
    mutados(1, :) = mejor_individuo;

    % Nueva generación
    poblacion = mutados;
end

% Mejor solución
final_aptitudes = evaluarAptitudHimmelblau(poblacion);
[~, idx_final] = max(final_aptitudes);
mejor = poblacion(idx_final, :);
costo_final = himmelblau(mejor(1), mejor(2));

fprintf('Mejor solución encontrada: x = %.4f, y = %.4f\n', mejor(1), mejor(2));
fprintf('Valor de la función: %.6f\n', costo_final);

% Gráfica de convergencia
figure;
plot(1:num_generaciones, maximos, 'LineWidth', 2);
xlabel('Generación');
ylabel('Valor mínimo de f(x, y)');
title('Convergencia - Himmelblau con GA');
grid on;

function aptitudes = evaluarAptitudHimmelblau(poblacion)
    [n, ~] = size(poblacion);
    aptitudes = zeros(n, 1);
    for i = 1:n
        x = poblacion(i,1);
        y = poblacion(i,2);
        f = himmelblau(x, y);
        aptitudes(i) = 1 / (1 + f); % maximizar
    end
end

function z = himmelblau(x, y)
    z = (x^2 + y - 11)^2 + (x + y^2 - 7)^2;
end

function seleccionados = seleccionPorTorneo(poblacion, fitness, k)
    len = size(poblacion, 1);
    seleccionados = zeros(size(poblacion));
    
    for i = 1:len
        indices = randsample(len, k);
        [~, best_idx] = max(fitness(indices));
        seleccionados(i,:) = poblacion(indices(best_idx), :);
    end
end

function hijos = cruceAritmetico(poblacion, prob)
    [n, d] = size(poblacion);
    hijos = zeros(n, d);
    for i = 1:2:n-1
        p1 = poblacion(i, :);
        p2 = poblacion(mod(i,n)+1, :);
        if rand < prob
            alpha = rand;
            hijos(i, :) = alpha*p1 + (1-alpha)*p2;
            hijos(i+1, :) = (1-alpha)*p1 + alpha*p2;
        else
            hijos(i, :) = p1;
            hijos(i+1, :) = p2;
        end
    end
end

function mutados = mutacionGaussiana(poblacion, prob, sigma)
    [n, d] = size(poblacion);
    ruido = sigma * randn(n, d);
    mascara = rand(n, d) < prob;
    mutados = poblacion;
    mutados(mascara) = mutados(mascara) + ruido(mascara);
end
