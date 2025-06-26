% Parámetros del problema (Antena)
area_size = 100;
r = 15; % radio de cobertura de cada antena
n_antenas = 5;

% Generar puntos de usuarios aleatoriamente
n_puntos = 300;
puntos = rand(n_puntos, 2) * area_size;

% Parámetros GA
tam_poblacion = 100;
num_generaciones = 200;
prob_cruce = 0.8;
prob_mutacion = 0.2;

% Cada individuo tiene 2 coordenadas por antena
num_genes = n_antenas * 2;
poblacion = rand(tam_poblacion, num_genes) * area_size;

mejores = zeros(num_generaciones, 1);

for gen = 1:num_generaciones
    aptitudes = evaluarCobertura(poblacion, puntos, r);

    % Guardar mejor
    [mejores(gen), idx] = max(aptitudes);
    elite = poblacion(idx,:);

    % Selección
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 5);

    % Cruce
    hijos = cruceDeUnPunto(seleccionados, prob_cruce);

    % Mutación
    mutados = mutacionGaussiana(hijos, prob_mutacion, 5, 0, area_size);

    % Elitismo
    mutados(1,:) = elite;
    poblacion = mutados;
end

% Mostrar mejor
[~, idx_final] = max(evaluarCobertura(poblacion, puntos, r));
mejor = poblacion(idx_final, :);
coords = reshape(mejor, 2, [])';

% Gráfica de cobertura
figure;
scatter(puntos(:,1), puntos(:,2), 'k.');
hold on;
viscircles(coords, r * ones(n_antenas,1), 'Color','b');
scatter(coords(:,1), coords(:,2), 100, 'r', 'filled');
title('Ubicación de Antenas y Cobertura');
axis([0 area_size 0 area_size]);
grid on;

% Convergencia
figure;
plot(1:num_generaciones, mejores, 'LineWidth', 2);
xlabel('Generación');
ylabel('Porcentaje de cobertura');
title('Convergencia del algoritmo genético');
grid on;

% ============ FUNCIONES AUXILIARES =============

function aptitudes = evaluarCobertura(poblacion, puntos, r)
    [n, d] = size(poblacion);
    n_antenas = d / 2;
    aptitudes = zeros(n,1);

    for i = 1:n
        coords = reshape(poblacion(i,:), 2, [])';
        cobertura = false(size(puntos,1), 1);

        for j = 1:n_antenas
            dist = sqrt(sum((puntos - coords(j,:)).^2, 2));
            cobertura = cobertura | (dist <= r);
        end

        aptitudes(i) = sum(cobertura) / size(puntos,1); % fracción cubierta
    end
end

function seleccionados = seleccionPorTorneo(poblacion, aptitudes, k)
    n = size(poblacion, 1);
    seleccionados = zeros(size(poblacion));
    for i = 1:n
        idx = randsample(n, k);
        [~, best] = max(aptitudes(idx));
        seleccionados(i,:) = poblacion(idx(best), :);
    end
end

function hijos = cruceDeUnPunto(poblacion, prob)
    [n, d] = size(poblacion);
    hijos = zeros(size(poblacion));
    for i = 1:2:n-1
        p1 = poblacion(i,:);
        p2 = poblacion(i+1,:);
        if rand < prob
            punto = randi([1 d-1]);
            hijos(i,:) = [p1(1:punto), p2(punto+1:end)];
            hijos(i+1,:) = [p2(1:punto), p1(punto+1:end)];
        else
            hijos(i,:) = p1;
            hijos(i+1,:) = p2;
        end
    end
end

function mutados = mutacionGaussiana(poblacion, prob, sigma, min_val, max_val)
    [n, d] = size(poblacion);
    mascara = rand(n, d) < prob;
    ruido = sigma * randn(n, d);
    mutados = poblacion;
    mutados(mascara) = mutados(mascara) + ruido(mascara);
    % Asegurar que se mantenga dentro del área
    mutados = max(mutados, min_val);
    mutados = min(mutados, max_val);
end

%% Parámetros del problema (Firewall)
num_atributos = 10; % número de condiciones binarias por paquete
num_reglas = 5;
long_regla = num_atributos;
long_individuo = num_reglas * long_regla;

% Dataset ficticio: 500 paquetes
n_paquetes = 500;
X = randi([0 1], n_paquetes, num_atributos);
Y = randi([0 1], n_paquetes, 1); % 1 = malicioso, 0 = benigno

% Parámetros del GA

tam_poblacion = 100;
num_generaciones = 200;
prob_cruce = 0.8;
prob_mutacion = 0.02;

poblacion = randi([0 1], tam_poblacion, long_individuo);
mejores = zeros(num_generaciones, 1);

for gen = 1:num_generaciones
    aptitudes = evaluarFirewall(poblacion, X, Y, num_reglas, long_regla);
    [mejores(gen), idx] = max(aptitudes);
    elite = poblacion(idx, :);

    % Selección por torneo
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 5);

    % Cruce de un punto
    hijos = cruceDeUnPunto(seleccionados, prob_cruce);

    % Mutación flip-bit
    mutados = mutacionFlipBit(hijos, prob_mutacion);

    % Elitismo
    mutados(1,:) = elite;
    poblacion = mutados;
end

% Mostrar resultados

[~, idx_final] = max(evaluarFirewall(poblacion, X, Y, num_reglas, long_regla));
mejor = poblacion(idx_final, :);
fprintf("Mejor individuo:\n%s\n", mat2str(mejor));
fprintf("Aptitud final: %.3f\n", mejores(end));

figure;
plot(1:num_generaciones, mejores, 'LineWidth', 2);
xlabel('Generación');
ylabel('Aptitud');
title('Optimización de reglas de firewall');
grid on;


% Funciones auxiliares


function aptitudes = evaluarFirewall(poblacion, X, Y, num_reglas, long_regla)
    [n, d] = size(poblacion);
    aptitudes = zeros(n, 1);
    n_paquetes = size(X,1);
    
    for i = 1:n
        reglas = reshape(poblacion(i,:), long_regla, [])';
        TP = 0; FP = 0;

        for j = 1:n_paquetes
            pkt = X(j,:);
            es_malicioso = Y(j);

            match = any(all(reglas == pkt, 2)); % alguna regla coincide

            if match && es_malicioso
                TP = TP + 1;
            elseif match && ~es_malicioso
                FP = FP + 1;
            end
        end

        aptitudes(i) = TP / (TP + FP + 1); % maximizamos
    end
end

function hijos = cruceDeUnPunto(poblacion, prob)
    [n, d] = size(poblacion);
    hijos = zeros(size(poblacion));
    for i = 1:2:n-1
        p1 = poblacion(i,:);
        p2 = poblacion(i+1,:);
        if rand < prob
            punto = randi([1 d-1]);
            hijos(i,:) = [p1(1:punto), p2(punto+1:end)];
            hijos(i+1,:) = [p2(1:punto), p1(punto+1:end)];
        else
            hijos(i,:) = p1;
            hijos(i+1,:) = p2;
        end
    end
end

function mutados = mutacionFlipBit(poblacion, prob)
    [n, d] = size(poblacion);
    mascara = rand(n, d) < prob;
    mutados = poblacion;
    mutados(mascara) = 1 - mutados(mascara);
end

function seleccionados = seleccionPorTorneo(poblacion, aptitudes, k)
    n = size(poblacion, 1);
    seleccionados = zeros(size(poblacion));
    for i = 1:n
        idx = randsample(n, k);
        [~, best] = max(aptitudes(idx));
        seleccionados(i,:) = poblacion(idx(best), :);
    end
end
