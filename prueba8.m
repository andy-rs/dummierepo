%% Parámetros del problema (Canales)
n_cells = 10;
n_channels = 15;

% Matriz de adyacencia
adyacencia = rand(n_cells) < 0.7;
adyacencia = triu(adyacencia, 1);
adyacencia = adyacencia + adyacencia';  % simétrica
adyacencia(1:n_cells+1:end) = 0;  % sin autointerferencia

% Parámetros del algoritmo genético
tam_poblacion = 100;
num_generaciones = 200;
prob_cruce = 0.8;
prob_mutacion = 0.1;

% Inicialización
poblacion = randi([1 n_channels], tam_poblacion, n_cells);
maximos = zeros(num_generaciones, 1);
mejores_individuos = zeros(num_generaciones, n_cells);

for gen = 1:num_generaciones
    % Evaluar aptitud (maximizar: menos interferencia)
    aptitudes = evaluarAptitudCanales(poblacion, adyacencia);
    interferencias = (1 ./ aptitudes) - 1;

    % Guardar mejor
    [~, idx] = max(aptitudes);
    mejor_individuo_actual = poblacion(idx, :);
    maximos(gen) = interferencias(idx);
    mejores_individuos(gen, :) = poblacion(idx, :);

    % Selección por torneo
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 10);

    % Cruce uniforme
    hijos = cruceUniforme(seleccionados, prob_cruce);

    % Mutación: cambio aleatorio de canal
    mutados = mutacionUniforme(hijos, prob_mutacion, 1, n_channels);

    % Elitismo
    mutados(1, :) = mejor_individuo_actual;

    % Siguiente generación
    poblacion = mutados;
end

% Mostrar mejor solución
final_aptitudes = evaluarAptitudCanales(poblacion, adyacencia);
[~, idx_final] = max(final_aptitudes);
mejor = poblacion(idx_final, :);
interferencia_final = contarInterferencias(mejor, adyacencia);
fprintf("Mejor asignación de canales: %s\n", mat2str(mejor));
fprintf("Interferencia total: %d\n", interferencia_final);

% Gráfica
figure;
plot(1:num_generaciones, maximos, 'LineWidth', 2);
xlabel("Generación");
ylabel("Interferencia del mejor individuo");
title("Evolución de interferencia en asignación de canales");
grid on;

function aptitudes = evaluarAptitudCanales(poblacion, adyacencia)
    [n, n_cells] = size(poblacion);
    aptitudes = zeros(n, 1);
    for i = 1:n
        individuo = poblacion(i, :);
        interferencias = contarInterferencias(individuo, adyacencia);
        aptitudes(i) = 1 / (1 + interferencias); % Maximizar
    end
end

function total = contarInterferencias(individuo, adyacencia)
    n = length(individuo);
    total = 0;
    for i = 1:n-1
        for j = i+1:n
            if adyacencia(i,j)
                diff = abs(individuo(i) - individuo(j));
                if diff == 0 || diff == 1
                    total = total + 1;
                end
            end
        end
    end
end

function mutados = mutacionUniforme(poblacion, prob, min_val, max_val)
    [len, num_genes] = size(poblacion);
    mascara = rand(len, num_genes) < prob;
    nuevos = randi([min_val max_val], len, num_genes);
    mutados = poblacion;
    mutados(mascara) = nuevos(mascara);
end

function hijos = cruceUniforme(poblacion, probabilidad)
    [len, num_genes] = size(poblacion);
    hijos = zeros(len, num_genes);
    for i = 1:2:len
        p1 = poblacion(i, :);
        p2 = poblacion(mod(i,len)+1, :);
        mask = rand(1, num_genes) < probabilidad;
        hijos(i,:) = p1;
        hijos(i+1,:) = p2;
        hijos(i,mask) = p2(mask);
        hijos(i+1,mask) = p1(mask);
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

%% Parámetros del escaneo (Puertos)
n_puertos = 20;
% GA
tam_poblacion = 200;
num_generaciones = 200;
prob_cruce = 0.8;
prob_mutacion = 0.15;
% Población inicial: permutaciones de puertos
poblacion = zeros(tam_poblacion, n_puertos);
for i = 1:tam_poblacion
    poblacion(i,:) = randperm(n_puertos);
end
mejores = zeros(num_generaciones,1);
for gen = 1:num_generaciones
    aptitudes = evaluarAptitudEscaneo(poblacion);
    [mejores(gen), idx] = max(aptitudes);
    elite = poblacion(idx,:);
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 5);
    hijos = cruceOX(seleccionados, prob_cruce);
    mutados = mutacionSwap(hijos, prob_mutacion);
    mutados(1,:) = elite; % elitismo
    poblacion = mutados;
end
% Mostrar resultado
[~, idx_final] = max(evaluarAptitudEscaneo(poblacion));
mejor = poblacion(idx_final, :);
fprintf("Mejor secuencia de escaneo:\n%s\n", mat2str(mejor));
figure;
plot(1:num_generaciones, mejores, 'LineWidth', 2);
xlabel('Generación');
ylabel('Aptitud (Inverso detección)');
title('Evasión en escaneo de puertos');
grid on;

% ----------------------
% Funciones auxiliares
% ----------------------

function aptitudes = evaluarAptitudEscaneo(poblacion)
    [n, m] = size(poblacion);
    aptitudes = zeros(n,1);
    for i = 1:n
        secuencia = poblacion(i,:);
        deteccion = scoreDeteccion(secuencia);
        aptitudes(i) = 1 / (1 + deteccion);
    end
end
function score = scoreDeteccion(secuencia)
    consecutivos = sum(abs(diff(secuencia)) == 1);
    cambios_grandes = sum(abs(diff(secuencia)) > 10);
    repetidos = length(unique(secuencia)) < length(secuencia);
    patron_secuencial = sum(diff(secuencia) > 0 & diff(secuencia) < 3);
    score = 3*consecutivos + 2*repetidos + 1*patron_secuencial - 0.3*cambios_grandes;
    score = max(0, score);
end
function hijos = cruceOX(poblacion, prob)
    [n, d] = size(poblacion);
    hijos = zeros(size(poblacion));
    for i = 1:2:n-1
        if i+1 <= n
            p1 = poblacion(i,:);
            p2 = poblacion(i+1,:);
            if rand < prob
                p1_cut = randi([2 d-2]);
                p2_cut = randi([p1_cut+1 d-1]);
                hijos(i,:) = ordenCruzado(p1, p2, p1_cut, p2_cut);
                hijos(i+1,:) = ordenCruzado(p2, p1, p1_cut, p2_cut);
            else
                hijos(i,:) = p1;
                hijos(i+1,:) = p2;
            end
        end
    end
end
function hijo = ordenCruzado(p1, p2, punto1, punto2)
    d = length(p1);
    hijo = zeros(1,d);
    hijo(punto1:punto2) = p1(punto1:punto2);
    resto = p2(~ismember(p2, hijo(punto1:punto2)));
    posiciones = [1:punto1-1, punto2+1:d];
    hijo(posiciones) = resto(1:length(posiciones));
end
function mutados = mutacionSwap(poblacion, prob)
    [n, d] = size(poblacion);
    mutados = poblacion;
    for i = 1:n
        if rand < prob
            idx = randperm(d, 2);
            mutados(i,[idx(1) idx(2)]) = mutados(i,[idx(2) idx(1)]);
        end
    end
end
function seleccionados = seleccionPorTorneo(poblacion, aptitudes, k)
    n = size(poblacion,1);
    seleccionados = zeros(size(poblacion));
    for i = 1:n
        idx = randi(n, k, 1);
        [~, best] = max(aptitudes(idx));
        seleccionados(i,:) = poblacion(idx(best), :);
    end
end
