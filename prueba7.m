% Parámetros del edificio (WIFI)
filas = 20; columnas = 20;
n_repetidores = 5;
radio = 4;

% Parámetros GA
tam_poblacion = 100;
num_generaciones = 200;
prob_cruce = 0.8;
prob_mutacion = 0.1;

% Longitud del cromosoma = 2 coordenadas por repetidor
long_individuo = n_repetidores * 2;

% Inicialización
poblacion = randi([1 max(filas,columnas)], tam_poblacion, long_individuo);
mejores = zeros(num_generaciones, 1);

for gen = 1:num_generaciones
    aptitudes = evaluarCobertura(poblacion, filas, columnas, radio, n_repetidores);
    [mejores(gen), idx] = max(aptitudes);
    elite = poblacion(idx, :);

    % Selección por torneo
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 5);

    % Cruce de un punto
    hijos = cruceDeUnPunto(seleccionados, prob_cruce);

    % Mutación uniforme
    mutados = mutacionUniforme(hijos, prob_mutacion, 1, max(filas,columnas));

    % Elitismo
    mutados(1,:) = elite;
    poblacion = mutados;
end

% Mostrar resultados
[~, idx_final] = max(evaluarCobertura(poblacion, filas, columnas, radio, n_repetidores));
mejor = poblacion(idx_final, :);
fprintf("Mejor ubicación de repetidores:\n");
disp(reshape(mejor, [2, n_repetidores])');

figure;
plot(1:num_generaciones, mejores, 'LineWidth', 2);
xlabel('Generación');
ylabel('Celdas cubiertas');
title('Cobertura Wi-Fi optimizada');
grid on;

% ----------------------------
% Funciones auxiliares
% ----------------------------

function aptitudes = evaluarCobertura(poblacion, filas, columnas, radio, n_repetidores)
    [n, d] = size(poblacion);
    aptitudes = zeros(n,1);

    for i = 1:n
        coords = reshape(poblacion(i,:), [2, n_repetidores])';
        grid = zeros(filas, columnas);
        for j = 1:n_repetidores
            x = round(coords(j,1));
            y = round(coords(j,2));
            if x < 1 || x > filas || y < 1 || y > columnas
                continue;
            end
            for dx = -radio:radio
                for dy = -radio:radio
                    if dx^2 + dy^2 <= radio^2
                        xi = x + dx;
                        yi = y + dy;
                        if xi >= 1 && xi <= filas && yi >= 1 && yi <= columnas
                            grid(xi, yi) = 1;
                        end
                    end
                end
            end
        end
        aptitudes(i) = sum(grid(:)); % total de celdas cubiertas
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

function mutados = mutacionUniforme(poblacion, prob, min_val, max_val)
    [n, d] = size(poblacion);
    mascara = rand(n, d) < prob;
    nuevos = randi([min_val max_val], n, d);
    mutados = poblacion;
    mutados(mascara) = nuevos(mascara);
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

%% Simulación: dataset sintético (Caracteristicas)
num_muestras = 1000;
num_features = 20;
X = rand(num_muestras, num_features);
y = randi([0 1], num_muestras, 1);  % clases binarias

% Parámetros GA
tam_poblacion = 50;
num_generaciones = 100;
prob_cruce = 0.8;
prob_mutacion = 0.05;

% Inicialización
poblacion = randi([0 1], tam_poblacion, num_features);
mejores = zeros(num_generaciones, 1);

for gen = 1:num_generaciones
    aptitudes = evaluarAptitudCiber(poblacion, X, y);

    % Guardar el mejor
    [mejores(gen), idx] = max(aptitudes);
    elite = poblacion(idx, :);

    % Selección
    seleccionados = seleccionPorTorneo(poblacion, aptitudes, 5);

    % Cruce
    hijos = cruceDeUnPunto(seleccionados, prob_cruce);

    % Mutación
    mutados = mutacionFlipBit(hijos, prob_mutacion);

    % Elitismo
    mutados(1,:) = elite;

    % Nueva generación
    poblacion = mutados;
end

% Mostrar resultado final
[~, idx_final] = max(evaluarAptitudCiber(poblacion, X, y));
mejor_individuo = poblacion(idx_final, :);
fprintf("Características seleccionadas: %s\n", mat2str(find(mejor_individuo)));
fprintf("Total: %d\n", sum(mejor_individuo));

% Gráfica
figure;
plot(1:num_generaciones, mejores, 'LineWidth', 2);
xlabel('Generación');
ylabel('Accuracy (fitness)');
title('Selección de características para detección de intrusos');
grid on;


% FUNCIONES AUXILIARES

function aptitudes = evaluarAptitudCiber(poblacion, X, y)
    [n, d] = size(poblacion);
    aptitudes = zeros(n,1);
    for i = 1:n
        mascara = poblacion(i,:) == 1;
        if sum(mascara) == 0
            aptitudes(i) = 0; % sin características => 0
        else
            X_sel = X(:,mascara);
            % Hold-out o CV simple (aquí, 70-30 hold-out con KNN)
            cv = cvpartition(size(X_sel,1), 'HoldOut', 0.3);
            modelo = fitcknn(X_sel(training(cv),:), y(training(cv)), 'NumNeighbors', 3);
            y_pred = predict(modelo, X_sel(test(cv),:));
            acc = mean(y_pred == y(test(cv)));
            aptitudes(i) = acc;
        end
    end
end

function seleccionados = seleccionPorTorneo(poblacion, aptitudes, k)
    n = size(poblacion,1);
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

function mutados = mutacionFlipBit(poblacion, prob)
    [n, d] = size(poblacion);
    mascara = rand(n, d) < prob;
    mutados = poblacion;
    mutados(mascara) = 1 - mutados(mascara);
end
