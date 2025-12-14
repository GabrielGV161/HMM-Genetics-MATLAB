%Queremos tres estados: Enlazante con el gen(1), No enlazante(2), Enlazante sin el gen(3)
%Probabilidades iniciales
PI = [0, 0.6, 0.1, 0.3]; %Si las distribuciones iniciales son iguales, al normalizar 
% no habrá diferencia entre EMIS y EMIS_ADJ
%Matrices de transición para cada caso
TRANS = [0.3, 0.3, 0.4; 
         0.1, 0.3, 0.6; 
         0.5, 0.2, 0.3]; 

TRANS_2 = [0.5, 0.2, 0.3; 
           0.2, 0.5, 0.3; 
           0.1, 0.4, 0.5];
%%
%Matrices de emisión
%Se hace que cada uno de los estados (filas) tenga una probabilidad de
%emisión en ciertos SNP dados los 20 que hay. Todas han de sumar 20, de ahí
%que se haga de 1, 10 a 1, 10 etc
EMIS = [0.1*ones(1, 10) 0.05*ones(1, 10);  % Estado 1 favorece los primeros 10 SNPs
        0.05*ones(1, 5) 0.1*ones(1, 15);   % Estado 2 favorece los últimos 15 SNPs
        0.05*ones(1, 20)];                 % Estado 3 emite todos los SNPs igual
%Los estados son independientes entre sí en términos de emisión, por eso
%las columnas no suman 1. Por ejemplo, el SNP nº1, tiene una probabilidad
%de 0.1 de emitirse en el estado 1, y un 0.9 de no emitirse en ese estado,
%todo independiente de su probabilidad de emisión en el estado 2.
%%
%Creamos matrices ampliadas para integrar PI
TRANS_HAT = [PI; zeros(size(TRANS,1),1) TRANS]; %Añado una fila con PI y una columna de ceros para que sea cuadrada
%Para TRANS_2
TRANS_HAT2 = [PI; zeros(size(TRANS_2,1),1) TRANS_2];
%Para EMIS
EMIS_HAT = [zeros(1,size(EMIS,2)); EMIS]; %Añado una fila de ceros para que sea compatible con TRANS_HAT


% Normalizar las matrices ampliadas
EMIS_HAT = EMIS_HAT ./ sum(EMIS_HAT, 2); % Normalizar la matriz EMIS_HAT
TRANS_HAT = TRANS_HAT ./ sum(TRANS_HAT, 2); % Normalizar la matriz TRANS_HAT
TRANS_HAT2 = TRANS_HAT2 ./ sum(TRANS_HAT2, 2); % Normalizar la matriz TRANS_HAT2
%Comprobación de que las matrices ampliadas son correctas
sum(TRANS_HAT, 2);  % Debe devolver 1 para todas las filas
sum(EMIS_HAT, 2); % Debe devolver 1 para todas las filas (excepto la fila inicial)
%%
%Gráficos de calor para visualizar las matrices
figure;
subplot(1,1,1);
heatmap(TRANS_HAT, 'Title', 'Matriz TRANS');
figure
subplot(1,1,1);
heatmap(EMIS_HAT, 'Title', 'Matriz de Emisión');
figure
subplot(1,1,1);
heatmap(TRANS_HAT2, 'Title', 'Matriz TRANS2');
figure
subplot(1,1,1);
heatmap(PI, 'Title', 'Vector de probabilidades iniciales')
%%
%Generamos secuencia de observaciones
n = input('Número de observaciones: '); %Longitud de la secuencia de observaciones (steps)
[seq, states] = hmmgenerate(n, TRANS_HAT, EMIS_HAT); %Generamos la secuencia
%La función hmmgenerate me va a generar una secuencia aleatoria de 1000 emisiones y estados
%teniendo en cuenta las matrices TRANS y EMIS.
%seq es la secuencia de 1000 observaciones (SNPS con valores del 1 al 20
%states es la secuencia de 1000 estados con valores del 2 al 4 porque el 1
%es el estado ficticio
%La función comienza en el estado 1 en el step 0 y transiciona al estado 1 en el step 0 por deafult
%%
%Calculamos la ruta más óptima con Viterbi
Vit1 = hmmviterbi(seq, TRANS_HAT, EMIS_HAT); %No usa PI porque asume que el modelo comienza en el estado 1
%Hay que multipilicar la matriz de EMIS por PI para quitar el sesgo del programa (ver más arriba)
Vit2 = hmmviterbi(seq, TRANS_HAT2, EMIS_HAT);
%%
%Calculamos con la foward variable la probabilidad de observar esta
%secuencia de observaciones generada. Recuerda que el foward modela la
%probabilidad de observar un símbolo en un estado t. De tal forma que si
%itera eso para todos los t, me da la probabilidad de una secuencia
%completa.
% En escala logarítmica dado que a secuencias muy grandes, la probabilidad tiende a 0
potito = pr_hmmMINE2(seq, TRANS_HAT, EMIS_HAT, PI); % Escala lineal con secuencias cortas
potito
% Calcular la matriz de probabilidades Forward
% Llamar a la función ajustada
[p, m] = pr_hmmMINE2(seq, TRANS_HAT, EMIS_HAT, PI);

% Excluir el estado ficticio
m_real = m(:, 2:end); % Eliminar la columna correspondiente al estado ficticio
%Gráfico de áreas apiladas de los distintos estados
%Definimos un alpha
[~, alpha] = pr_hmmMINE2(seq, TRANS_HAT, EMIS_HAT, PI);
alpha_normalized = alpha ./ sum(alpha, 2); %Lo normalizamos
%Se puede identificar regiones donde un estado domina analizando alpha a lo
%largo de los SNPs. Recuerda que alpha en foward me da una probabilidad
%para cada estado.
figure;
area(alpha_normalized(:, 2:end), 'LineStyle', 'none'); % Excluye el estado ficticio
xlabel('SNP');
ylabel('Probabilidad normalizada');
title('Distribución de probabilidades por estado');
legend({'Enlazante con el gen', 'No enlazante', 'Enlazante sin el gen'}, 'Location', 'Best');
colormap('parula');
grid on;
%Gráfico de líneas
figure;
hold on;
plot(alpha_normalized(:, 2), '-o', 'LineWidth', 2, 'DisplayName', 'Enlazante con el gen');
plot(alpha_normalized(:, 3), '-s', 'LineWidth', 2, 'DisplayName', 'No enlazante');
plot(alpha_normalized(:, 4), '-^', 'LineWidth', 2, 'DisplayName', 'Enlazante sin el gen');
xlabel('SNP');
ylabel('Probabilidad normalizada');
title('Evolución de probabilidades por estado');
legend show;
grid on;
hold off;
%Gráfico de calor compacto
figure;
imagesc(alpha_normalized(:, 2:end)');
colormap('hot');
colorbar;
xlabel('SNP');
ylabel('Estados');
yticks(1:3);
yticklabels({'Enlazante con el gen', 'No enlazante', 'Enlazante sin el gen'});
title('Probabilidades normalizadas por estado');
%%
%Lo mismo para TRANS2
% Calcular la matriz de probabilidades Forward
% Llamar a la función ajustada
[p, m] = pr_hmmMINE2(seq, TRANS_HAT2, EMIS_HAT, PI);

% Excluir el estado ficticio
m_real = m(:, 2:end); % Eliminar la columna correspondiente al estado ficticio
%Gráfico de áreas apiladas de los distintos estados
%Definimos un alpha
[~, alpha2] = pr_hmmMINE2(seq, TRANS_HAT2, EMIS_HAT, PI);
alpha_normalized = alpha2 ./ sum(alpha2, 2); %Lo normalizamos
%Se puede identificar regiones donde un estado domina analizando alpha a lo
%largo de los SNPs. Recuerda que alpha en foward me da una probabilidad
%para cada estado.
figure;
area(alpha_normalized(:, 2:end), 'LineStyle', 'none'); % Excluye el estado ficticio
xlabel('SNP');
ylabel('Probabilidad normalizada');
title('Distribución de probabilidades por estado TRANS2');
legend({'Enlazante con el gen', 'No enlazante', 'Enlazante sin el gen'}, 'Location', 'Best');
colormap('parula');
grid on;
%Gráfico de líneas
figure;
hold on;
plot(alpha_normalized(:, 2), '-o', 'LineWidth', 2, 'DisplayName', 'Enlazante con el gen');
plot(alpha_normalized(:, 3), '-s', 'LineWidth', 2, 'DisplayName', 'No enlazante');
plot(alpha_normalized(:, 4), '-^', 'LineWidth', 2, 'DisplayName', 'Enlazante sin el gen');
xlabel('SNP');
ylabel('Probabilidad normalizada');
title('Evolución de probabilidades por estado TRANS2');
legend show;
grid on;
hold off;
%Gráfico de calor compacto
figure;
imagesc(alpha_normalized(:, 2:end)');
colormap('hot');
colorbar;
xlabel('SNP');
ylabel('Estados');
yticks(1:3);
yticklabels({'Enlazante con el gen', 'No enlazante', 'Enlazante sin el gen'});
title('Probabilidades normalizadas por estado TRANS2');
%%
%Ploteamos Vit1 con diagrama de calor
%
figure;
subplot(1,1,1);
heatmap(Vit1, 'Title', 'Matriz de Vit1');
%Vit2
figure;
subplot(1,1,1);
heatmap(Vit2, 'Title', 'Matriz de Vit2');
%Hay que buscar una manera de comparar Viterbi con otras secuencias
%plotteando la probabilidida de cada step en el diagrama
%%
%Diagrama para la secuencia de estados Vit1
% Crear matriz para el diagrama de calor
n_steps = length(Vit1); % Longitud de la secuencia
heatmap_data = zeros(3, n_steps); % Tres estados reales (excluye el estado ficticio)

for t = 1:n_steps
    heatmap_data(Vit1(t)-1, t) = 1; % Resta 1 para mapear estados 2, 3, 4 a filas 1, 2, 3
end

% Crear el gráfico de calor con imágenes
figure;
imagesc(heatmap_data);

% Cambiar el esquema de colores
colormap('jet'); % Esquema de colores más vistoso
colorbar; % Agrega una barra de colores para referencia

% Ajustar los ejes
xlabel('Step');
ylabel('Estado');
title('Ruta de estados (Viterbi1)');

% Ajustar las etiquetas del eje Y
yticks(1:3); % Fila 1 a 3
yticklabels({'Enlazante con el gen', 'No enlazante', 'Enlazante sin el gen'});
ax = gca; % Obtener el objeto del eje
ax.YAxis.TickLabelInterpreter = 'none'; % Asegurar que se respeten las etiquetas

%Diagrama para la secuencia de estados Vit2
% Crear matriz para el diagrama de calor
n_steps = length(Vit2); % Longitud de la secuencia
heatmap_data = zeros(3, n_steps); % Tres estados reales (excluye el estado ficticio)

for t = 1:n_steps
    heatmap_data(Vit2(t)-1, t) = 1; % Resta 1 para mapear estados 2, 3, 4 a filas 1, 2, 3
end

% Crear el gráfico de calor con imágenes
figure;
imagesc(heatmap_data);

% Cambiar el esquema de colores
colormap('jet'); % Esquema de colores más vistoso
colorbar; % Agrega una barra de colores para referencia

% Ajustar los ejes
xlabel('Step');
ylabel('Estado');
title('Ruta de estados (Viterbi2)');

% Ajustar las etiquetas del eje Y
yticks(1:3); % Fila 1 a 3
yticklabels({'Enlazante con el gen', 'No enlazante', 'Enlazante sin el gen'});
ax = gca; % Obtener el objeto del eje
ax.YAxis.TickLabelInterpreter = 'none'; % Asegurar que se respeten las etiquetas
%%
%Comparación de rutas
% Generar posiciones ficticias
positions = 1:length(Vit1); % Suponemos que ambas rutas tienen la misma longitud

% Diagrama de calor para comparar rutas
comparison_data = [Vit1; Vit2];
figure;
imagesc(comparison_data);
colormap('parula');
colorbar;
xlabel('Posición en el genoma');
ylabel('Configuraciones');
yticks(1:2);
yticklabels({'Conf TRANS', 'Conf TRANS2'});
title('Comparación de rutas (diagrama de calor)');
