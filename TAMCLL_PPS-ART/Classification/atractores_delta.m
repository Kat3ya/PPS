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
for i1 = 1:length(De)
    Attr(i1,:) = dec2bin(De(i1),26);
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
%%
i = 2^(26);
for m = 1:length(De)
    j(:,m) =(i-De(m));
end
j=j';
%%
Fenotipo = cellstr(X{:,27});
J = array2table(j);
tabla = [J, X(:,27)];
M0 = find(strcmp(Fenotipo,'M0'));
M1 = find(strcmp(Fenotipo,'M1'));
M2 = find(strcmp(Fenotipo,'M2'));
NLC = find(strcmp(Fenotipo,'NLC'));
%%
for n1 = 1:length(M0)
    M0_delta(:,n1) = j(M0(n1));
end
M0_delta = array2table(M0_delta');

for n2 = 1:length(M1)
    M1_delta(:,n2) = j(M1(n2));
end
M1_delta = array2table(M1_delta');

for n3 = 1:length(M2)
    M2_delta(:,n3) = j(M2(n3));
end
M2_delta = array2table(M2_delta');

for n4 = 1:length(NLC)
    NLC_delta(:,n4) = j(NLC(n4));
end
NLC_delta = array2table(NLC_delta');
%% De, clase
tabla2 = [A1(:,4) X(:,27)];