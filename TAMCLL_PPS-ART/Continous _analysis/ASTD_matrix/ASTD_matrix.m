close all
clear all
clc
%% Load all files
% Archivo con lista de atractores y clases (fenotipos)
Attr = load("TAMCLL_attractors_class.mat");
Attractors = Attr.X;
%Archivo .csv con transiciones entre atrctores, cambia entre OR o AND
Attr2 = readtable("Tabla_de_transiciones_1-26_AND.csv");       
%Archivo .txt con lista de genes que generan transiciones, cambia entre OR o AND
Genes2 = textscan(fopen('headers2_AND.txt'),'%s %s','Delimiter',',');       
Genes2 = Genes2{1};
%Ciclo for para obtener la matriz ASTD (Available State Transitions Distribution) de cada gen 
for i=1:length(Genes2)
    Genes = Attr2{:,56};
    z = find(strcmp(Genes, Genes2{i}));         %Para encontrar cuántas veces se repite ese gen
    z1 = length(z);
    xi = Attr2{z(1):z(end),2:27};
    xf = Attr2{z(1):z(end),29:54};

    % Tiempo de manipulación promedio
    Phenotypes = Attr2{z(1):z(end),1};
    %M2
    zM2 = find(strcmp(Phenotypes, 'M2'));
    if isempty(zM2)
        tprom_M2 = 0;
    else
        tprom_M2 = mean(Attr2{zM2(1):zM2(end),55});
    end
    %M1
    zM1 = find(strcmp(Phenotypes, 'M1'));
    if isempty(zM1)
        tprom_M1 = 0;
    else
        tprom_M1 = mean(Attr2{zM1(1):zM1(end),55});
    end
    %NLC
    zNLC = find(strcmp(Phenotypes, 'NLC')); 
    if isempty(zNLC)
        tprom_NLC = 0;
    else
        tprom_NLC = mean(Attr2{zNLC(1):zNLC(end),55});
    end
    %M0
    zM0 = find(strcmp(Phenotypes, 'M0')); 
    if isempty(zM0)
        tprom_M0 = 0;
    else
        tprom_M0 = mean(Attr2{zM0(1):zM0(end),55});
    end
    %% Matriz de transiciones 1384x1384
    % Convertimos a equivalente decimal para facilitar la búsqueda de valores
    % Lista de 1384 atractores
    De = table2array(Attractors(:,1:26));
    for i1 = 1:length(De)
        Attr_de(i1,:) = bin2dec(num2str(De(i1,:)));
    end
    % Lista de equivalentes decimales de estados iniciales y finales que sí transicionan
    for i2 = 1:length(xi)
        xi_de(i2,:) = bin2dec(num2str(xi(i2,:)));
        xf_de(i2,:) = bin2dec(num2str(xf(i2,:)));
    end
    % En esta parte se genera una matriz (Transistions_m) de zeros de 1384x1384 para luego
    % reemplazar con 1 en los estados que haya transición
    Transitions_m = zeros(1384,1384);                   
    for i3 = 1:length(xi_de)
        x_0 = find(Attr_de==xi_de(i3));
        x_d = find(Attr_de==xf_de(i3));
        Transitions_m(x_0,x_d) = 1;
    end
    % Se obtiene la matriz donde 0-no transiciona, 1-sí transiciona
    %figure
    %heatmap(Transitions_m,'Xlabel','Initial State x_0','Ylabel','Final State x_d','Title',['Transitions of TAMCLL-GRN via ',Genes2{i},' (AND operator)'],'ColorbarVisible','off','Colormap',hot)
    fprintf('El gen %s genera %d transiciones con operación OR\n',Genes2{i},z1);
    fprintf('El tiempo de manipulación promedio para M2 es %g\n',tprom_M2);
    fprintf('El tiempo de manipulación promedio para M1 es %g\n',tprom_M1);
    fprintf('El tiempo de manipulación promedio para NLC es %g\n',tprom_NLC);
    fprintf('El tiempo de manipulación promedio para M0 es %g\n',tprom_M0);
    fprintf('------------------------------------------------\n');
    %% Matriz ASTD
    % Seccionamos la matriz Transitions_m para obtener las concurrencias de
    % la transición entre fenotipos 
    Cprop = PropChart([320 78 60 926]',Transitions_m);
    xvalues = {'M2','M1','NLC','M0'};
    yvalues = {'M2','M1','NLC','M0'};
    figure
    heatmap(xvalues,yvalues,Cprop,'Xlabel','Initial State x_0','Ylabel','Final State x_d','Title',['Transitions of TAMCLL-GRN via ',Genes2{i},' (AND operator)'],'ColorLimits',[0 100])
end





