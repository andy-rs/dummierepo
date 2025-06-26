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

% Cruce OX
function hijos = cruceOX(poblacion)
    len = size(poblacion, 1);
    num_genes = size(poblacion, 2);
    hijos = zeros(len, num_genes);

    for i = 1:2:len
        p1 = poblacion(i, :);
        p2 = poblacion(i+1, :);
        puntos = sort(randsample(num_genes, 2));

        c1 = zeros(1, num_genes);
        c2 = zeros(1, num_genes);
        c1(puntos(1):puntos(2)) = p1(puntos(1):puntos(2));
        c2(puntos(1):puntos(2)) = p2(puntos(1):puntos(2));

        idx1 = puntos(2) + 1;
        idx2 = puntos(2) + 1;
        fill_pos1 = mod(puntos(2), num_genes) + 1;
        fill_pos2 = fill_pos1;

        for k = 1:num_genes
            g2 = p2(mod(idx1-1, num_genes) + 1);
            if ~ismember(g2, c1)
                c1(fill_pos1) = g2;
                fill_pos1 = mod(fill_pos1, num_genes) + 1;
            end
            idx1 = idx1 + 1;

            g1 = p1(mod(idx2-1, num_genes) + 1);
            if ~ismember(g1, c2)
                c2(fill_pos2) = g1;
                fill_pos2 = mod(fill_pos2, num_genes) + 1;
            end
            idx2 = idx2 + 1;
        end

        hijos(i, :) = c1;
        hijos(i+1, :) = c2;
    end
end
