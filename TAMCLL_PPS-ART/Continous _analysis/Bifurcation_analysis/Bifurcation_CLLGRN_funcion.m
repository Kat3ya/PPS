% This code computes the bifurcation diagrams of the neccesary time (k) of the
% time window of each controller to produce the transition between
% attractors
clear all
clc
%% Load all files
Attr = load("TAMCLL_attractors_class.mat");
Attractors = Attr.X;
De = table2array(Attractors(:,1:26));
for i = 1:length(De)
    Attr_de(i,:) = bin2dec(num2str(De(i,:)));
end
%% Boolean equations TAMCLL-GRN
% IL1(t+1) = IL1
% IL10(t+1) = IL10
% HMGB1(t+1) = HMGB1
% IFNg(t+1) = IFNg
% GMCSF(t+1) = GMCSF
% LPS(t+1) = LPS
% IC(t+1) = IC
% IL13(t+1) = IL13
% IL4(t+1) = IL4
% MCSF(t+1) = MCSF
% IFNgR(t+1) = IFNg | (LPS & ! (IC & (LPS | IL1))) & ! (STAT6 | STAT1)
% STAT1(t+1) = IFNgR | STAT1 & ! STAT6
% STAT5(t+1) = GMCSF & ! (STAT3 | IRF4)
% NFkb(t+1) = (STAT1 | TNFa | (LPS & ! (IC & (LPS | IL1))) | (IL1 | (NFkb | TNFa))) 
%   & ! (STAT6 | (IC & (LPS | IL1)) | PPARg | STAT6)
% PPARg(t+1) = (IL4 & IL13) | MCSF | (IC & (LPS | IL1)) & ! STAT6
% STAT6(t+1) = (IL4 & IL13) | MCSF
% JMJD3(t+1) = (IL4 & IL13) | MCSF
% STAT3(t+1) = ((IL10 | ((PPARg | STAT3) & ! (IRF5 | TNFa))) | EGF | STAT3) & ! ((IC & (LPS | IL1)) | PPARg)
% IRF4(t+1) = JMJD3
% IRF5(t+1) = STAT5 & ! IRF4
% EGF(t+1) = (IC & (LPS | IL1)) | STAT3
% IL12(t+1) = STAT1 | STAT5 | NFkb
% TNFa(t+1) = IRF5 & !((PPARg | STAT3) & ! (IRF5 | TNFa))
% TGFb(t+1) = STAT3 & ! TNFa
% HIF1A(t+1) = (STAT3 | ((PPARg | STAT3) & ! (IRF5 | TNFa))) & ! STAT1
% RAGE(t+1) = HMGB1
%% Attractors
% Initial condition 1384 attractors
i_c = Attractors{:,1:26};
for m = 1:size(i_c,1)
    x_c(m)=sum(i_c(m,:));  % Sum of active/inactive nodes for the initial condition (Value between 0-26)
end
%% Iterative variation of parameter k
Ci = [];
for l=1:size(i_c,1)
    [X_i,X_f,x_c1] = Bifurcation_CLL(i_c(l,:),x_c(l));
    Xi(l,:) = X_i(1,:);
    Xf(l,:) = X_f(end,:);
    Xc(:,l) = x_c1;
end
%%
contador = 0;
for n=1:size(i_c,1)
    G1 = bin2dec(num2str(Xi(n,:)));
    G2 = bin2dec(num2str(Xf(n,:)));
    if G1 ~= G2
        B=[];
        contador = contador + 1;
        ci = Xi(n,:);
        cf = Xf(n,:);
        xc = Xc(:,n);

        Ci(contador,:) = ci;
        Cf(contador,:) = cf;
        Ci_De(contador) = bin2dec(num2str(Ci(contador,:)));
        class1(contador) = find(Attr_de == Ci_De(contador));
        Fenotipo1(contador,:) = Attractors(class1(contador),27);
        Cf_De(contador) = bin2dec(num2str(Cf(contador,:)));
        class2(contador) = find(Attr_de == Cf_De(contador));
        Fenotipo2(contador,:) = Attractors(class2(contador),27);
        %disp('sí transiciona')

        X_c(:,contador) = round(xc);
        B = find(X_c(:,contador)==X_c(end,contador));
        tt(contador,:) = B(1)*0.01;
        
    else
        %disp('no transiciona')
    end
    
end
Fenotipo1 = table2array(Fenotipo1);
Fenotipo2 = table2array(Fenotipo2);

Data = table(Fenotipo1, Ci, Fenotipo2, Cf, tt);
Data = renamevars(Data,"Ci","Estado inicial");
Data = renamevars(Data,"Cf","Estado final");
Data = renamevars(Data,"tt","t de manipulación");

Data2 = table(Fenotipo1, Fenotipo2, tt);

