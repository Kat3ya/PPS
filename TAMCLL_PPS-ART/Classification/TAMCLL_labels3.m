clc
clear all
%% 
% load all files
Genes = textscan(fopen('headers2.txt'),'%s %s','Delimiter',',');
A1 = readtable('TAM1384.csv');
load('svm3_model.mat')
load('kNN3_model.mat')
%% 
De = table2array(A1(:,4));
for i = 1:length(De)
    Attr(i,:) = dec2bin(De(i),26);
end
AttractorsTable = array2table(Attr);
AttractorsTable.Properties.VariableNames = Genes{1,1};
%% Assign labels to data to categorize them
sortedNames = sort(AttractorsTable.Properties.VariableNames(1:26)); % Reorder table columns in alphabetical order
X = AttractorsTable(:,sortedNames);
T = X{:,1:26}; % Create a matrix of data for manipulation
% Assign phenotype to each attractor according to study
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
%%
Z = [A2(:,2) A2(:,5)];
yfit = kNN3_model.predictFcn(Z);
class = categorical(X{:,27});
h = confusionchart(class,yfit);
h.Title = {'kNN model'};
A3 = [AttractorsTable X(:,27) table(yfit)];
A3 = renamevars(A3,"Var27","Labels");
A3 = renamevars(A3,"yfit","Predicted");
%% Supervised Classification
Z = [A2(:,2) A2(:,5)];
yfit = kNN3_model.predictFcn(Z);
%yfit = smv3_model.predictFcn(Z);
class = categorical(X{:,27});
h = confusionchart(class,yfit);
h.Title = {'kNN model'};
A3 = [AttractorsTable X(:,27) table(yfit)];
A3 = renamevars(A3,"Var27","Labels");
A3 = renamevars(A3,"yfit","Predicted");
for j=1:1384
   if class(j) == yfit(j)
       A4(j,:) = zeros(1,26);
   elseif class(j) ~= yfit(j)
       A4(j,:) = T(j,1:26);
   end
end
suma = sum(A4,2);
f = find(suma);
for l = 1:length(f)
    A5(l,:)=A4(f(l),:);
    yfit2(l) = yfit(f(l));
    class2(l) = class(f(l));
    B(l,:) = A2(f(l),5); %Interaccioes M2 predichas
end
yfit2=yfit2';
class2=class2';
%%
A6 = array2table(A5);
A6.Properties.VariableNames = Genes{1,1};
etiquetas=table(class2, yfit2);
A6 = [A6 etiquetas];
T2 = A5(:,1:26);
[~,s]=sortrows(A6(:,28));
A6 = A6(s,:);
%%
%{
A7 = movevars(A6,'IL12','Before',1);
A8 = movevars(A7,'NFkb','After',1);
A9 = movevars(A8,'TNFa','After',2);
A10 = movevars(A9,'STAT1','After',3);
A11 = movevars(A10,'STAT5','After',4);
A12 = movevars(A11,'IL10','After',5);
A13 = movevars(A12,'STAT3','After',6);
A14 = movevars(A13,'STAT6','After',7);
A15 = movevars(A14,'PPARg','After',8);
A16 = movevars(A15,'TGFb','After',9);
A17 = movevars(A16,'EGF','After',10);
A18 = movevars(A17,'RAGE','After',11);
%}
%%
genes = A6.Properties.VariableNames;
h1 = figure(1);
j1 = heatmap(1:37,genes(1:26),T2');
j1.Title = "Expression profiles of recovered attractors";
j1.ColorbarVisible = 'off';
j1.GridVisible = 'off';
print(h1,'Attractors','-dmeta','-r1000');
%% Perfil de expresión clases 
cat = categories(A6.class2);
K = [categorical(0); cat];
y = categorical(AttractorsTable.Properties.VariableNames);
A = zeros(4,26);
oj= [];
for t = 2:5
    i = A6.class2 == K(t);
    r = A6(i,1:26);
    u = zeros(height(r),26);
    for k = 1:height(r)
        u(k,:) = r{k,:};
    end
    A(t-1,:) = mean(u);
    oj= [oj, A(t-1,:)];
    figure
    bar(y,A(t-1,:))
end
%% Perfil de expresión clases predichas
cat = categories(A6.yfit2);
K = [categorical(0); cat];
y = categorical(AttractorsTable.Properties.VariableNames);
A = zeros(4,26);
oj= [];
for t = 2:5
    i = A6.yfit2 == K(t);
    r = A6(i,1:26);
    u = zeros(height(r),26);
    for k = 1:height(r)
        u(k,:) = r{k,:};
    end
    A(t-1,:) = mean(u);
    oj= [oj, A(t-1,:)];
    figure
    bar(y,A(t-1,:))
end
%% Unsupervised classification
% k means attractor, basin size
Y1 = table2array(A2(:,2:7));
Y = [Y1(:,1) Y1(:,2)];
[cluster_k,C] = kmeans(Y,4);
figure
gscatter(Y(:,1),Y(:,2),cluster_k);
title('K-Means Using Squared Euclidean Distance')
hold on
plot(C(:,1),C(:,2),'kx')
xlabel('Attractor')
ylabel('Basin size')

labels = table2array(A2(:,8));
cluster_k = categorical(cluster_k);
figure
c1= confusionmat(labels,cluster_k);
confusionchart(c1)
title('K-Means Clustering with Squared Euclidean Distance')

%{
figure
plot(Y(cluster_k==1,1),Y(cluster_k==1,2),'r.','MarkerSize',12)
hold on
plot(Y(cluster_k==2,1),Y(cluster_k==2,2),'b.','MarkerSize',12)
plot(C(:,1),C(:,2),'kx','MarkerSize',15,'LineWidth',3) 
legend('Cluster 1','Cluster 2','Centroids')
title 'Cluster Assignments and Centroids'
xlabel('Basin size')
ylabel('De')
hold off
%}

% k means attractor, De
Y2 = [Y1(:,1) Y1(:,3)];
[cluster_k2,C2] = kmeans(Y2,4);
figure
gscatter(Y2(:,1),Y2(:,2),cluster_k2);
title('K-Means Using Squared Euclidean Distance')
hold on
plot(C2(:,1),C2(:,2),'kx')
xlabel('Attractor')
ylabel('De')

cluster_k2 = categorical(cluster_k2);
figure
c2= confusionmat(labels,cluster_k2);
confusionchart(c2)
title('K-Means Clustering with Squared Euclidean Distance')

% k means attractor, interactions
Y3 = [Y1(:,1) Y1(:,4)];
[cluster_k3,C3] = kmeans(Y3,4);
figure
gscatter(Y3(:,1),Y3(:,2),cluster_k3);
title('K-Means Using Squared Euclidean Distance')
hold on
plot(C3(:,1),C3(:,2),'kx')
xlabel('Attractor')
ylabel('Interactions')

cluster_k3 = categorical(cluster_k3);
figure
c3= confusionmat(labels,cluster_k3);
confusionchart(c3)
title('K-Means Clustering with Squared Euclidean Distance')

% k means attractor, c. largo
Y4 = [Y1(:,1) Y1(:,5)];
[cluster_k4,C4] = kmeans(Y4,4);
figure
gscatter(Y4(:,1),Y4(:,2),cluster_k4);
title('K-Means Using Squared Euclidean Distance')
hold on
plot(C4(:,1),C4(:,2),'kx')
xlabel('Attractor')
ylabel('C.largo')

cluster_k4 = categorical(cluster_k4);
figure
c4= confusionmat(labels,cluster_k4);
confusionchart(c4)
title('K-Means Clustering with Squared Euclidean Distance')

% k means hamming distance
Y5 = table2array(A2(:,2));
[cluster_k5,C5] = kmeans(Y5,4);
%{
figure
gscatter(Y5(:,1),Y5(:,2),cluster_k5);
title('K-Means Using Squared Euclidean Distance')
hold on
plot(C5(:,1),C5(:,2),'kx')
xlabel('Attractor')
ylabel('C.largo')
%}

cluster_k5 = categorical(cluster_k5);
figure
c5= confusionmat(labels,cluster_k5);
confusionchart(c5)
title('K-Means Clustering with Hamming Distance')
%%
Y6 = table2array(A1(:,2:6));
Y6 = [Y6(:,3) Y6(:,4)];

[cluster_k6,C6] = kmeans(Y6,4);

figure
gscatter(Y6(:,1),Y6(:,2),cluster_k6);
title('K-Means Using Squared Euclidean Distance')
hold on
plot(C6(:,1),C6(:,2),'kx')
xlabel('Basin size')
ylabel('Interactions')


labels = double(table2array(X(:,27)));
cluster_k6 = categorical(cluster_k6);
figure
c6 = confusionmat(labels,cluster_k6);
confusionchart(c6)
title('K-Means Clustering')
%% Histograma modelos
%SVM
%modelo = categorical({'SVM1', 'kNN1', 'SVM2', 'kNN2', 'SVM3', 'kNN3', 'SVM4', 'kNN4'});
%acc = [0 68.1; 68.1 88.4; 68.1 72.5; 68.1 93.5; 68.1 81.2];
%acc = [68.1 20.3; 67.0 21; 68.1 4.4; 67.0 0; 68.1 25.4; 67.0 28.7; 68.1 13.1; 67.0 9.4];
modelo = categorical({'1', '2', '3', '4'});
acc = [88.3 88; 72.5 67.0; 93.5 95.7; 81.2 76.4];
g = bar(modelo,acc);
xtips1 = g(1).XEndPoints;
ytips1 = g(1).YEndPoints;
labels1 = string(g(1).YData);
text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

xtips2 = g(2).XEndPoints;
ytips2 = g(2).YEndPoints;
labels2 = string(g(2).YData);
text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xlabel('Modelo')
ylabel('Porcentaje exactitud (%)')
title('Comparación modelos aprendizaje supervisado')

%%
%kNN
figure
modelo2 = categorical({'kNN1', 'kNN2', 'kNN3', 'kNN4'});
acc2 = [67.0 21; 67.0 0; 67.0 28.7; 67.0 9.4];
v = bar(modelo2,acc2,'stacked');

xtips2v = v(2).XEndPoints;
ytips2v = v(2).YEndPoints;
labels2v = string(v(2).YData);
text(xtips2v,ytips2v,labels2v,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
xlabel('Modelo')
ylabel('Porcentaje exactitud (%)')
title('Comparación modelos kNN')
%%
C=[A6(20,:);A6(26,:);A6(27,:);A6(30,:);A6(33,:);A6(34,:);A6(37,:)];
cat = categories(C.yfit2);
K = [categorical(0); categorical(2)];
y = categorical(AttractorsTable.Properties.VariableNames);
A = zeros(4,26);
oj= [];
for t = 2:5
    i = C.yfit2 == K(t);
    r = C(i,1:26);
    u = zeros(height(r),26);
    for k = 1:height(r)
        u(k,:) = r{k,:};
    end
    A(t-1,:) = mean(u);
    oj= [oj, A(t-1,:)];
    figure
    bar(y,A(t-1,:))
end
C2=[B(20,:);B(26,:);B(27,:);B(30,:);B(33,:);B(34,:);B(37,:)];
mean(C2);
%%
[~,c]=sortrows(A2(:,8));
B1 = A2(c,:);
%%
figure
gscatter(Z(:,1),Z(:,2),cat);
%% Set de entrenamiento
% Comprobamos que el algoritmo toma el 80% de cada clase para el set de 
% entrenamiento total del 80% de los 1384 atractores
h = kNN3_model.ClassificationKNN.Y;
h = double(string(h));
M0_kNN = find(h==4);
M1_kNN = find(h==1);
M2_kNN = find(h==2);
NLC_kNN = find(h==3);

g = smv3_model.ClassificationSVM.Y;
g = double(string(g));
M0_SVM = find(h==4);
M1_SVM = find(h==1);
M2_SVM = find(h==2);
NLC_SVM = find(h==3)