% Parámetros del problema (Reinas)
N = 8; 

% Parámetros del algoritmo genético
num_generaciones = 200;
tam_poblacion = 100;
prob_cruce = 0.8;
prob_mutacion = 0.2;

% Inicialización
poblacion = randi([1 N], tam_poblacion, N);
minimos = zeros(num_generaciones, 1);
mejores_individuos = zeros(num_generaciones, N);

for gen = 1:num_generaciones
    % Evaluar aptitud
    aptitudes = evaluarAptitudReinas(poblacion);
    conflictos = (1 ./ aptitudes) - 1;

    % Guardar mejor individuo
    [~, idx] = max(aptitudes);
    minimos(gen) = conflictos(idx)
    mejores_individuos(gen, :) = poblacion(idx, :);

    % Selección por ranking
    seleccionados = seleccionPorRanking(poblacion, aptitudes);

    % Cruce uniforme
    hijos = cruceUniforme(seleccionados, prob_cruce);

    % Mutación por fila aleatoria entera
    mutados = mutacionUniforme(hijos, prob_mutacion, 1, N);

    % Elitismo
    mutados(1, :) = poblacion(idx, :);

    % Nueva generación
    poblacion = mutados;
end

% Mostrar mejor solución encontrada
aptitudes_finales = evaluarAptitudReinas(poblacion);
[~, idx_mejor] = max(aptitudes_finales);
mejor_individuo = poblacion(idx_mejor, :);
conflictos_final = contarConflictos(mejor_individuo);
fprintf('Mejor individuo: %s\n', mat2str(mejor_individuo));
fprintf('Conflictos: %d\n', conflictos_final);

% Gráfica de convergencia
figure;
plot(1:num_generaciones, minimos, 'LineWidth', 2);
xlabel('Generación');
ylabel('Número mínimo de conflictos');
title(sprintf('Evolución de conflictos para %d reinas', N));
grid on;

% ========= Funciones =========

function aptitudes = evaluarAptitudReinas(poblacion)
    [n, N] = size(poblacion);
    aptitudes = zeros(n, 1);
    for i = 1:n
        conflictos = contarConflictos(poblacion(i, :));
        aptitudes(i) = 1 / (1 + conflictos); 
    end
end

function total = contarConflictos(individuo)
    N = length(individuo);
    total = 0;
    for i = 1:N-1
        for j = i+1:N
            if individuo(i) == individuo(j)
                total = total + 1;
            elseif abs(individuo(i) - individuo(j)) == abs(i - j)
                total = total + 1;
            end
        end
    end
end

function seleccionados = seleccionPorRanking(poblacion, fitness)
    len = size(poblacion, 1);
    [~, orden] = sort(fitness, 'descend');
    rangos = len:-1:1;
    probs = rangos / sum(rangos);
    indices = randsample(orden, len, true, probs);
    seleccionados = poblacion(indices, :); 
end

function hijos = cruceUniforme(poblacion, probabilidad)
    [len, num_genes] = size(poblacion);
    hijos = zeros(len, num_genes);
    for i = 1:2:len
        p1 = poblacion(i, :);
        p2 = poblacion(mod(i,len)+1, :);
        mask = rand(1, num_genes) < probabilidad;
        hijos(i, :) = p1;
        hijos(i+1, :) = p2;
        hijos(i, mask) = p2(mask);
        hijos(i+1, mask) = p1(mask);
    end
end

function mutados = mutacionUniforme(poblacion, prob, min_val, max_val)
    [len, num_genes] = size(poblacion);
    mascara = rand(len, num_genes) < prob;
    nuevos = randi([min_val max_val], len, num_genes);
    mutados = poblacion;
    mutados(mascara) = nuevos(mascara);
end

%% Datos a ajustar (Parabola)
x = linspace(-10, 10, 100);
y_real = 2*x.^2 + 3*x + 5 + randn(size(x))*10;

% Parámetros del problema
n_coef = 3;  % a, b, c
rango_min = -10;
rango_max = 10;

% Parámetros GA
tam_poblacion = 100;
num_generaciones = 200;
prob_cruce = 0.8;
prob_mutacion = 0.1;

% Inicialización
poblacion = rango_min + (rango_max - rango_min) * rand(tam_poblacion, n_coef);
maximos = zeros(num_generaciones, 1);
mejores_individuos = zeros(num_generaciones, n_coef);

for gen = 1:num_generaciones
    % Evaluar aptitud
    aptitudes = evaluarAptitudParabola(poblacion, x, y_real);
    errores = (1 ./ aptitudes) - 1;

    % Guardar mejor
    [~, idx] = max(aptitudes);
    maximos(gen) = errores(idx);
    mejores_individuos(gen, :) = poblacion(idx, :);
    mejor_individuo_actual = poblacion(idx, :);

    % Selección por torneo
    seleccionados = seleccionPorRuleta(poblacion, aptitudes);

    % Cruce aritmético
    hijos = cruceAritmetico(seleccionados, prob_cruce);

    % Mutación gaussiana
    mutados = mutacionGaussiana(hijos, prob_mutacion, 1);

    % Elitismo
    mutados(1, :) = mejor_individuo_actual;

    % Nueva generación
    poblacion = mutados;
end

% Mostrar mejor resultado
final_aptitudes = evaluarAptitudParabola(poblacion, x, y_real);
[~, idx_final] = max(final_aptitudes);
mejor = poblacion(idx_final, :);
fprintf("Mejor ajuste encontrado: a = %.4f, b = %.4f, c = %.4f\n", mejor);
y_pred = mejor(1)*x.^2 + mejor(2)*x + mejor(3);

% Gráfica de resultados
figure;
plot(x, y_real, 'b.', 'DisplayName', 'Datos');
hold on;
plot(x, y_pred, 'r-', 'LineWidth', 2, 'DisplayName', 'Ajuste GA');
legend;
xlabel('x');
ylabel('y');
title('Ajuste de parábola con algoritmo genético');
grid on;

% Gráfica de convergencia
figure;
plot(1:num_generaciones, maximos, 'LineWidth', 2);
xlabel('Generación');
ylabel('Error mínimo (ECM)');
title('Convergencia del algoritmo');
grid on;

% ============== FUNCIONES ==============

function aptitudes = evaluarAptitudParabola(poblacion, x, y)
    [n, ~] = size(poblacion);
    aptitudes = zeros(n, 1);
    for i = 1:n
        a = poblacion(i,1);
        b = poblacion(i,2);
        c = poblacion(i,3);
        y_pred = a*x.^2 + b*x + c;
        ecm = mean((y - y_pred).^2);
        aptitudes(i) = 1 / (1 + ecm);  % maximización
    end
end

function seleccionados = seleccionPorRuleta(poblacion, fitness)
    len = size(poblacion, 1);
    probs = fitness / sum(fitness);
    indices_seleccionados = randsample(1:len, len, true, probs);
    seleccionados = poblacion(indices_seleccionados, :);
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

