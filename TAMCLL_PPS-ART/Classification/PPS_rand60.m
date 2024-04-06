clc
clear all
%% CLASIFICACIÓN DE ATRACTORES DE TAMCLL BGRN 
% En este código se seleccionan aleatoriamente 60 atractores de cada clase 
% y se entrenan modelos de aprendizaje supervisado para clasificarlos en 
% M0, M1, M2 y NLC.
%% Importamos el archivo de atractores y carcaterísticas
Genes = textscan(fopen('headers2.txt'),'%s %s','Delimiter',',');
A1 = readtable('TAM1384.csv');
load('rand5_60_kNNmodel.mat')
load('rand7_60_kNNmodel.mat')
load('kNN3_model.mat')
%% Obtenemos el equivalente decimal y hacemos una tabla con todas las caracteristícas
De = table2array(A1(:,4));
for i = 1:length(De)
    Attr(i,:) = dec2bin(De(i),26);
end
AttractorsTable = array2table(Attr);
AttractorsTable.Properties.VariableNames = Genes{1,1};
%% Asignamos las etiquetas a los datos pra categorizarlos Assign 
sortedNames = sort(AttractorsTable.Properties.VariableNames(1:26)); % Se reordenan las columnas de la tabla en orden alfabetico
X = AttractorsTable(:,sortedNames);
T = X{:,1:26}; % Se crea una matriz de los datos para su manipulación
% Asignamos cada atractor a su fenotipo de acuerdo al artículo
for k = 1:1384
    if (T(k,10) == 1) && (T(k,18) == 1) && (T(k,26) == 1) && ((T(k,21) == 1) || (T(k,23) == 1))
        X{k,27} = categorical("M1");
    elseif (T(k,9) == 1) && ((T(k,22) == 1) || (T(k,24) == 1)) && (T(k,19) == 1)
        X{k,27} = categorical("M2");
    elseif (T(k,25) == 1) && (T(k,3) == 1) && (T(k,1) == 1) && (T(k,20) == 1)
        X{k,27} = categorical("NLC");
    else
        X{k,27} = categorical("M0");
    end
end
% 926 atractores en M0 66%
% 78 atractores en M1 5.6%
% 320 atractores en M2 23.1%
% 60 atractores en NLC 4.3%
LOGattr = log10(table2array(A1(:,2)));
A2 = [A1 array2table(LOGattr) X(:,27)];
A2 = renamevars(A2,"Var27","Labels");
%A2(1,:) = [];
%% Seleccionamos aleatoriamente 60 atractores de cada clase
% Se reordena la tabla por clase en orden 4, 3, 2, 1
[~,s]=sortrows(A2(:,8));
A3 = A2(s,:);
% Se separan los atractores por clase
class_M0 = A3(1:926,:);
class_NLC = A3(927:986,:);
class_M1 = A3(987:1064,:);
class_M2 = A3(1065:1384,:);
% Selección aleatoria de atactores
r_M0 = datasample(class_M0,60);
r_M1 = datasample(class_M1,60);
r_M2 = datasample(class_M2,60);
%% Modelo de clasificación con aprendizaje supervisado
% Se crea la tabla de datos para implentar el modelo SVM y KNN
A4 = [r_M0; r_M1; r_M2; class_NLC];
%% Clasificación 1384 atractores
Z = [A2(:,2) A2(:,5)];
Z1 = kNN3_model.ClassificationKNN.X;
%yfit = rand5_60_kNNmodel.predictFcn(Z1);
yfit = rand7_60_kNNmodel.predictFcn(Z1);
%yfit = kNN3_model.predictFcn(Z);
%class = categorical(X{:,27});
class = kNN3_model.ClassificationKNN.Y;
h = confusionchart(class,yfit);
h.Title = {'kNN model'};




