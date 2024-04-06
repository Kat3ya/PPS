clc
clear all 
close all
%% Boolean dynamic simulation TAMCLL GRN
% This code computes de dynamic simulation of a Boolean network based on
% its Boolean functions
%% Description of network
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
%% Initial conditions
%x0 = [1 0 0 0 1 0 0 0 0 0 1 1 0 0 0 1 0 0 1 1 0 0 0 0 0 0];    %Initial condition to M1, De = 35701952
%x0 = [0 1 1 1 0 1 1 1 1 1 0 0 0 0 1 1 1 0 1 0 1 0 0 0 1 1];    %Initial condition to M2, De = 31395491
%x0 = [0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 1 1 1];    %Initial condition to NLC, De = 25166375
x0 = [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0];    %Initial condition to M0, De = 16777280
IL1 = x0(:,1);
IL10 = x0(:,2);
HMGB1 = x0(:,3);
IFNg = x0(:,4);
GMCSF = x0(:,5);
LPS = x0(:,6);
IC = x0(:,7);
IL13 = x0(:,8);
IL4 = x0(:,9);
MCSF = x0(:,10);
IFNgR = x0(:,11);
STAT1 = x0(:,12);
STAT5 = x0(:,13);
NFkb = x0(:,14);
PPARg = x0(:,15);
STAT6 = x0(:,16);
JMJD3 = x0(:,17);
STAT3 = x0(:,18);
IRF4 = x0(:,19);
IRF5 = x0(:,20);
EGF = x0(:,21);
IL12 = x0(:,22);
TNFa = x0(:,23);
TGFb = x0(:,24);
HIF1A = x0(:,25);
RAGE = x0(:,26);
%% Parameters of simulation
i = 1;                  % Index of simulation
s = 7;                 % number of step
IL1(:,i) = IL1;
IL10(:,i) = IL10;
HMGB1(:,i) = HMGB1;
IFNg(:,i) = IFNg;
GMCSF(:,i) = GMCSF;
LPS(:,i) = LPS;
IC(:,i) = IC;
IL13(:,i) = IL13;
IL4(:,i) = IL4;
MCSF(:,i) = MCSF;
IFNgR(:,i) = IFNgR;
STAT1(:,i) = STAT1;
STAT5(:,i) = STAT5;
NFkb(:,i) = NFkb;
PPARg(:,i) = PPARg;
STAT6(:,i) = STAT6;
JMJD3(:,i) = JMJD3;
STAT3(:,i) = STAT3;
IRF4(:,i) = IRF4;
IRF5(:,i) = IRF5;
EGF(:,i) = EGF;
IL12(:,i) = IL12;
TNFa(:,i) = TNFa;
TGFb(:,i) = TGFb;
HIF1A(:,i) = HIF1A;
RAGE(:,i) = RAGE;
xs = [IL1(:,1), IL10(:,1), HMGB1(:,1), IFNg(:,1), GMCSF(:,1), LPS(:,1),...
    IC(:,1), IL13(:,1), IL4(:,1), MCSF(:,1), IFNgR(:,1), STAT1(:,1),...
    STAT5(:,1), NFkb(:,1), PPARg(:,1), STAT6(:,1), JMJD3(:,1),...
    STAT3(:,1), IRF4(:,1), IRF5(:,1), EGF(:,1), IL12(:,1), TNFa(:,1),...
    TGFb(:,1), HIF1A(:,1), RAGE(:,1)];
%% Boolean simulation of the network
for i=1:s
     IL1(:,i+1) = IL1(:,i);
     IL10(:,i+1) = IL10(:,i);
     HMGB1(:,i+1) = HMGB1(:,i);
     IFNg(:,i+1) = IFNg(:,i);
     GMCSF(:,i+1) = GMCSF(:,i);
     LPS(:,i+1) = LPS(:,i);
     IC(:,i+1) = IC(:,i);
     IL13(:,i+1) = IL13(:,i);
     IL4(:,i+1) = IL4(:,i);
     MCSF(:,i+1) = MCSF(:,i);
     IFNgR(:,i+1) = IFNg(:,i) | (LPS(:,i) & not(IC(:,i) & (LPS(:,i) |...
         IL1(:,i)))) & not(STAT6(:,i) | STAT1(:,i));
     STAT1(:,i+1) = IFNgR(:,i) | STAT1(:,i) & not(STAT6(:,i));
     STAT5(:,i+1) = GMCSF(:,i) & not(STAT3(:,i) | IRF4(:,i));
     NFkb(:,i+1) = (STAT1(:,i) | TNFa(:,i) | (LPS(:,i) & not(IC(:,i) &...
         (LPS(:,i) | IL1(:,i)))) | (IL1(:,i) | (NFkb(:,i) | TNFa(:,i))))...
         & not(STAT6(:,i) | (IC(:,i) & (LPS(:,i) | IL1(:,i))) | PPARg(:,i) | STAT6(:,i));
     PPARg(:,i+1) = (IL4(:,i) & IL13(:,i)) | MCSF(:,i) | (IC(:,i) & (LPS(:,i) | IL1(:,i))) & not(STAT6(:,i));
     STAT6(:,i+1) = (IL4(:,i) & IL13(:,i)) | MCSF(:,i);
     JMJD3(:,i+1) = (IL4(:,i) & IL13(:,i)) | MCSF(:,i);
     STAT3(:,i+1) = ((IL10(:,i) | ((PPARg(:,i) | STAT3(:,i)) & not(IRF5(:,i)...
         | TNFa(:,i)))) | EGF(:,i) | STAT3(:,i)) & not((IC(:,i) & (LPS(:,i)...
         | IL1(:,i))) | PPARg(:,i));
     IRF4(:,i+1) = JMJD3(:,i);
     IRF5(:,i+1) = STAT5(:,i) & not(IRF4(:,i));
     EGF(:,i+1) = (IC(:,i) & (LPS(:,i) | IL1(:,i))) | STAT3(:,i);
     IL12(:,i+1) = STAT1(:,i) | STAT5(:,i) | NFkb(:,i);
     TNFa(:,i+1) = IRF5(:,i) & not((PPARg(:,i) | STAT3(:,i)) & not(IRF5(:,i) | TNFa(:,i)));
     TGFb(:,i+1) = STAT3(:,i) & not(TNFa(:,i));
     HIF1A(:,i+1) = (STAT3(:,i) | ((PPARg(:,i) | STAT3(:,i)) & not(IRF5(:,i) | TNFa(:,i)))) & not(STAT1(:,i));
     RAGE(:,i+1) = HMGB1(:,i);
     xt(i,:) = [IL1(:,i+1), IL10(:,i+1), HMGB1(:,i+1), IFNg(:,i+1),...
         GMCSF(:,i+1), LPS(:,i+1), IC(:,i+1), IL13(:,i+1), IL4(:,i+1),...
         MCSF(:,i+1), IFNgR(:,i+1), STAT1(:,i+1), STAT5(:,i+1),...
         NFkb(:,i+1), PPARg(:,i+1), STAT6(:,i+1), JMJD3(:,i+1), STAT3(:,i+1),...
         IRF4(:,i+1), IRF5(:,i+1), EGF(:,i+1), IL12(:,i+1), TNFa(:,i+1),...
         TGFb(:,i+1), HIF1A(:,i+1), RAGE(:,i+1)];
end
%% Visualization of dyynamic transition of the network
xt
% Binario a decimal
r = bin2dec(num2str([x0;xt]))
% Decimal a binario 
% dec2bin(6,3)
stairs(r)
