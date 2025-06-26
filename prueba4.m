% Parámetros del problema (Job Sequencing)
n_jobs = 10;
profits = [100 19 27 25 15 90 30 10 40 70];
deadlines = [2 1 2 1 3 5 7 9 3 2];

% Parámetros GA
tam_poblacion = 200;
num_generaciones = 500;
prob_cruce = 0.7;
prob_mutacion = 0.1;

% Inicialización: población de permutaciones
poblacion = zeros(tam_poblacion, n_jobs);
for i = 1:tam_poblacion
    poblacion(i,:) = randperm(n_jobs);
end

maximos = zeros(num_generaciones, 1);
mejores_individuos = zeros(num_generaciones, n_jobs);

for gen = 1:num_generaciones
    % Evaluar aptitud (ganancia total)
    aptitudes = evaluarAptitudJobs(poblacion, profits, deadlines);
    [~, idx] = max(aptitudes);
    maximos(gen) = aptitudes(idx);
    mejores_individuos(gen,:) = poblacion(idx,:);
    mejor_individuo = poblacion(idx,:);

    % Selección por torneo
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 3);

    % Cruce ordenado (OX)
    hijos = cruceOX(seleccionados, prob_cruce);

    % Mutación por intercambio
    mutados = mutacionSwap(hijos, prob_mutacion);

    % Elitismo
    mutados(1,:) = mejor_individuo;

    % Siguiente generación
    poblacion = mutados;
end

% Mostrar mejor solución
final_aptitudes = evaluarAptitudJobs(poblacion, profits, deadlines);
[~, idx_final] = max(final_aptitudes);
mejor = poblacion(idx_final, :);
fprintf('Mejor orden de trabajos: %s\n', mat2str(mejor));
fprintf('Ganancia total: %d\n', evaluarAptitudJobs(mejor, profits, deadlines));

% Gráfica de convergencia
figure;
plot(1:num_generaciones, maximos, 'LineWidth', 2);
xlabel('Generación');
ylabel('Ganancia total');
title('Convergencia - Job Sequencing Problem');
grid on;

function aptitudes = evaluarAptitudJobs(poblacion, profits, deadlines)
    [n, n_jobs] = size(poblacion);
    aptitudes = zeros(n, 1);

    for i = 1:n
        orden = poblacion(i,:);
        ocupado = false(1, max(deadlines));
        ganancia = 0;

        for j = 1:n_jobs
            job = orden(j);
            d = deadlines(job);
            for t = d:-1:1
                if ~ocupado(t)
                    ocupado(t) = true;
                    ganancia = ganancia + profits(job);
                    break;
                end
            end
        end

        aptitudes(i) = ganancia; % ya es maximización
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

function hijos = cruceOX(poblacion, prob)
    [n, d] = size(poblacion);
    hijos = zeros(size(poblacion));
    for i = 1:2:n-1
        p1 = poblacion(i,:);
        p2 = poblacion(i+1,:);
        if rand < prob
            punto1 = randi([1 d-1]);
            punto2 = randi([punto1+1 d]);
            h1 = ordenCruzado(p1, p2, punto1, punto2);
            h2 = ordenCruzado(p2, p1, punto1, punto2);
            hijos(i,:) = h1;
            hijos(i+1,:) = h2;
        else
            hijos(i,:) = p1;
            hijos(i+1,:) = p2;
        end
    end
end

function hijo = ordenCruzado(p1, p2, punto1, punto2)
    d = length(p1);
    hijo = zeros(1, d);
    hijo(punto1:punto2) = p1(punto1:punto2);
    pos = mod(punto2, d) + 1;
    j = pos;
    for k = 1:d
        gene = p2(mod(punto2 + k - 1, d) + 1);
        if ~ismember(gene, hijo)
            hijo(j) = gene;
            j = mod(j, d) + 1;
        end
    end
end

function mutados = mutacionSwap(poblacion, prob)
    [n, d] = size(poblacion);
    mutados = poblacion;
    for i = 1:n
        if rand < prob
            idx = randperm(d, 2);
            tmp = mutados(i,idx(1));
            mutados(i,idx(1)) = mutados(i,idx(2));
            mutados(i,idx(2)) = tmp;
        end
    end
end

%% Parámetros del problema (TSP)
n_ciudades = 20;
coordenadas = rand(n_ciudades, 2) * 100;

% Parámetros del GA
tam_poblacion = 100;
num_generaciones = 200;
prob_cruce = 0.8;
prob_mutacion = 0.2;

% Inicialización de la población (permutaciones)
poblacion = zeros(tam_poblacion, n_ciudades);
for i = 1:tam_poblacion
    poblacion(i,:) = randperm(n_ciudades);
end

maximos = zeros(num_generaciones, 1);
mejores_individuos = zeros(num_generaciones, n_ciudades);

for gen = 1:num_generaciones
    aptitudes = evaluarAptitudTSP(poblacion, coordenadas);

    [~, idx] = max(aptitudes);
    maximos(gen) = 1 / aptitudes(idx) - 1;
    mejores_individuos(gen, :) = poblacion(idx,:);
    mejor_individuo = poblacion(idx,:);

    % Selección
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 5);

    % Cruce ordenado
    hijos = cruceOX(seleccionados, prob_cruce);

    % Mutación swap
    mutados = mutacionSwap(hijos, prob_mutacion);

    % Elitismo
    mutados(1,:) = mejor_individuo;

    % Nueva generación
    poblacion = mutados;
end

% Mejor ruta final
final_aptitudes = evaluarAptitudTSP(poblacion, coordenadas);
[~, idx_final] = max(final_aptitudes);
mejor = poblacion(idx_final, :);
mejor_ruta = [mejor mejor(1)];  % vuelta a la ciudad inicial
fprintf("Distancia final: %.2f\n", 1 / final_aptitudes(idx_final) - 1);

% Gráfica
figure;
plot(1:num_generaciones, maximos, 'LineWidth', 2);
xlabel('Generación');
ylabel('Distancia mínima');
title('Convergencia - TSP');
grid on;

% Ruta final
figure;
plot(coordenadas(mejor_ruta,1), coordenadas(mejor_ruta,2), '-o');
title('Mejor Ruta Encontrada');
xlabel('X');
ylabel('Y');
axis equal;
grid on;


function aptitudes = evaluarAptitudTSP(poblacion, coords)
    [n, d] = size(poblacion);
    aptitudes = zeros(n, 1);
    for i = 1:n
        ruta = poblacion(i,:);
        ruta = [ruta ruta(1)];
        dist = 0;
        for j = 1:d
            a = coords(ruta(j), :);
            b = coords(ruta(j+1), :);
            dist = dist + norm(a - b);
        end
        aptitudes(i) = 1 / (1 + dist); % maximización
    end
end

function seleccionados = seleccionPorTorneo(poblacion, aptitudes, k)
    n = size(poblacion, 1);
    seleccionados = zeros(size(poblacion));
    for i = 1:n
        idxs = randsample(n, k);
        [~, best] = max(aptitudes(idxs));
        seleccionados(i,:) = poblacion(idxs(best), :);
    end
end

function hijos = cruceOX(poblacion, prob)
    [n, d] = size(poblacion);
    hijos = zeros(size(poblacion));
    for i = 1:2:n-1
        p1 = poblacion(i,:);
        p2 = poblacion(i+1,:);
        if rand < prob
            punto1 = randi([1 d-1]);
            punto2 = randi([punto1+1 d]);
            h1 = ordenCruzado(p1, p2, punto1, punto2);
            h2 = ordenCruzado(p2, p1, punto1, punto2);
            hijos(i,:) = h1;
            hijos(i+1,:) = h2;
        else
            hijos(i,:) = p1;
            hijos(i+1,:) = p2;
        end
    end
end

function hijo = ordenCruzado(p1, p2, punto1, punto2)
    d = length(p1);
    hijo = zeros(1, d);
    hijo(punto1:punto2) = p1(punto1:punto2);
    pos = mod(punto2, d) + 1;
    j = pos;
    for k = 1:d
        gene = p2(mod(punto2 + k - 1, d) + 1);
        if ~ismember(gene, hijo)
            hijo(j) = gene;
            j = mod(j, d) + 1;
        end
    end
end

function mutados = mutacionSwap(poblacion, prob)
    [n, d] = size(poblacion);
    mutados = poblacion;
    for i = 1:n
        if rand < prob
            idx = randperm(d, 2);
            temp = mutados(i, idx(1));
            mutados(i, idx(1)) = mutados(i, idx(2));
            mutados(i, idx(2)) = temp;
        end
    end
end
