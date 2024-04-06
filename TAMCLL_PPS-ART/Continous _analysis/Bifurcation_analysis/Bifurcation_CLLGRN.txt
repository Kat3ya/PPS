% This code computes the bifurcation diagrams of the neccesary time (k) of the
% time window of each controller to produce the transition between
% attractors
clear all
clc
%% Load all files
A = textscan(fopen('headers2.txt'),'%s %s','Delimiter',',');
Genes = A{1}';
Attr = load("TAMCLL_attractors_class.mat");
Attractors = Attr.X;
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
%{
i_c = Attractors{:,1:26};
for m = 1:length(i_c)
    x_c(m,:)=sum(i_c(m));
end
%}
i_c=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %Initial condition
x_c=sum(i_c); % Sum of active/inactive nodes for the initial condition (Value between 0-3)
%% Iterative variation of parameter k
for k=0:0.01:3
    i1=1;
    for tu=0:0.01:10
        if (tu>0 && tu<=k)
             u(:,i1)=1;  %High level
             i1=i1+1;
        else
            u(:,i1)=0;  %Low level
            i1=i1+1;
        end
    end
    tu=0:0.01:10; %Active time of controller
    %% Solving ODEs
    options = odeset('RelTol',1e-3,'AbsTol',[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4...
        1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4...
        1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4]); %Configuration of ODE45 solver
    [T,X] = ode45(@(t,x) B2Diff_CLLGRN(t,x,tu,u),[0 10],i_c,options); %ODE45 solver 
    x_1=X(:,1); x_2=X(:,2); x_3=X(:,3); x_4=X(:,4); x_5=X(:,5); x_6=X(:,6);
    x_7=X(:,7); x_8=X(:,8); x_9=X(:,9); x_10=X(:,10); x_11=X(:,11); x_12=X(:,12);
    x_13=X(:,13); x_14=X(:,14); x_15=X(:,15); x_16=X(:,16); x_17=X(:,17); x_18=X(:,18);
    x_19=X(:,19); x_20=X(:,20); x_21=X(:,21); x_22=X(:,22); x_23=X(:,23); x_24=X(:,24);
    x_25=X(:,25); x_26=X(:,26);
    X_i=[x_1 x_2 x_3 x_4 x_5 x_6 x_7 x_8 x_9 x_10 x_11 x_12 x_13 x_14 x_15 x_16 x_17 x_18 x_19 x_20 x_21 x_22 x_23 x_24...
        x_25 x_26];
    %Sum of active/inactive nodes after solvig differential equation (Value between 0-13)
    x_s=x_1(end)+x_2(end)+x_3(end)+x_4(end)+x_5(end)+x_6(end)+x_7(end)+x_8(end)+x_9(end)+x_10(end)+x_11(end)+x_12(end)...
        +x_13(end)+x_14(end)+x_15(end)+x_16(end)+x_17(end)+x_18(end)+x_19(end)+x_20(end)+x_21(end)+x_22(end)+x_23(end)...
        +x_24(end)+x_25(end)+x_26(end); 
    x_c=[x_c,x_s];
    
    %% Visualization of solved differential equation system
    
    subplot(7,4,1);
    plot(T,X(:,1),'o-.')
    title('IL1')
    ylim([0 1]);
    subplot(7,4,2);
    plot(T,X(:,2),'o-.')
    title('IL10')
    ylim([0 1]);
    subplot(7,4,3);
    plot(T,X(:,3),'o-.')
    title('HMGB1')
    ylim([0 1]);
    subplot(7,4,4);
    plot(T,X(:,4),'o-.')
    title('IFNg')
    ylim([0 1]);
    subplot(7,4,5);
    plot(T,X(:,5),'o-.')
    title('GMCSF')
    ylim([0 1]);
    subplot(7,4,6);
    plot(T,X(:,6),'o-.')
    title('LPS')
    ylim([0 1]);
      subplot(7,4,7);
    plot(T,X(:,7),'o-.')
    title('IC')
    ylim([0 1]);
    subplot(7,4,8);
    plot(T,X(:,8),'o-.')
    title('IL13')
    ylim([0 1]);
    subplot(7,4,9);
    plot(T,X(:,9),'o-.')
    title('IL4')
    ylim([0 1]);
    subplot(7,4,10);
    plot(T,X(:,10),'o-.')
    title('MCSF')
    ylim([0 1]);
    subplot(7,4,11);
    plot(T,X(:,11),'o-.')
    title('IFNgR')
    ylim([0 1]);
    subplot(7,4,12);
    plot(T,X(:,12),'o-.')
    title('STAT1')
    ylim([0 1]);
    subplot(7,4,13);
    plot(T,X(:,13),'o-.')
    title('STAT5')
    ylim([0 1]);
    subplot(7,4,14);
    plot(T,X(:,14),'o-.')
    title('NFkb')
    ylim([0 1]);
    subplot(7,4,15);
    plot(T,X(:,15),'o-.')
    title('PPARg')
    ylim([0 1]);
    subplot(7,4,16);
    plot(T,X(:,16),'o-.')
    title('STAT6')
    ylim([0 1]);
    subplot(7,4,17);
    plot(T,X(:,17),'o-.')
    title('JMJD3')
    ylim([0 1]);
    subplot(7,4,18);
    plot(T,X(:,18),'o-.')
    title('STAT3')
    ylim([0 1]);
    subplot(7,4,19);
    plot(T,X(:,19),'o-.')
    title('IRF4')
    ylim([0 1]);
    subplot(7,4,20);
    plot(T,X(:,20),'o-.')
    title('IRF5')
    ylim([0 1]);
    subplot(7,4,21);
    plot(T,X(:,21),'o-.')
    title('EGF')
    ylim([0 1]);
    subplot(7,4,22);
    plot(T,X(:,22),'o-.')
    title('IL12')
    ylim([0 1]);
    subplot(7,4,23);
    plot(T,X(:,23),'o-.')
    title('TNFa')
    ylim([0 1]);
    subplot(7,4,24);
    plot(T,X(:,24),'o-.')
    title('TGFb')
    ylim([0 1]);
    subplot(7,4,25);
    plot(T,X(:,25),'o-.')
    title('HIF1A')
    ylim([0 1]);
    subplot(7,4,26);
    plot(T,X(:,26),'o-.')
    title('RAGE')
    ylim([0 1]);
    pause(0.1);
end
figure 
plot(x_c','o')
X_f=[round(x_1(end)) round(x_2(end)) round(x_3(end)) round(x_4(end)) round(x_5(end)) round(x_6(end)) round(x_7(end))...
    round(x_8(end)) round(x_9(end)) round(x_10(end)) round(x_11(end)) round(x_12(end)) round(x_13(end)) round(x_14(end))...
    round(x_15(end)) round(x_16(end)) round(x_17(end)) round(x_18(end)) round(x_19(end)) round(x_20(end)) round(x_21(end))...
    round(x_22(end)) round(x_23(end)) round(x_24(end)) round(x_25(end)) round(x_26(end))];