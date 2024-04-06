clc
clear all 
close all
%% Boolean dynamic simulation TAMCLL GRN
% This code computes de dynamic simulation of a Boolean network based on
% its Boolean functions
%% Description of network
% IFNg(t+1)   = IFNg
% GMCSF(t+1)  = GMCSF
% IL1(t+1)    = IL1
% LPS(t+1)    = LPS
% IC(t+1)     = IC
% IL13(t+1)   = IL13
% IL4(t+1)    = IL4
% IL10(t+1)   = IL10
% MCSF(t+1)   = MCSF
% HMGB1(t+1)  = HMGB1
% IFNgR(t+1)  = IFNg | IFNab & ! SOCS1
% CSF2Ra(t+1) = GMCSF
% IL1R(t+1)   = IL1 | IL1b
% TLR4(t+1)   = LPS & ! FCgR
% FCgR(t+1)   = IC & (LPS | IL1)
% IL4Ra(t+1)  = IL4 & IL13
% IL10r(t+1)  = IL10 | IL10s
% MCSFr(t+1)  = MCSF
% STAT1(t+1)  = IFNgR | STAT1 & ! STAT6
% STAT5(t+1)  = CSF2Ra & ! (STAT3 | IRF4)
% NFkb(t+1)   = (STAT1 | TNFa | TLR4 | IL1R) & ! (STAT6 | FCgR | PPARg | KLF4)
% PPARg(t+1)  = IL4Ra | MCSFr | ERK & ! STAT6
% STAT6(t+1)  = IL4Ra | MCSFr
% JMJD3(t+1)  = IL4Ra | MCSFr
% STAT3(t+1)  = (IL10r | EGF | STAT3) & ! (FCgR | PPARg)
% IRF3(t+1)   = TLR4
% ERK(t+1)    = FCgR
% KLF4(t+1)   = STAT6
% SOCS1(t+1)  = STAT6 | STAT1
% IRF4(t+1)   = JMJD3
% IRF5(t+1)   = STAT5 & ! IRF4
% IL1b(t+1)   = NFkb | TNFa
% IFNab(t+1)  = IRF3
% EGF(t+1)    = ERK | STAT3
% IL12(t+1)   = STAT1 | STAT5 | NFkb
% IL10s(t+1)  = (PPARg | STAT3) & ! (IRF5 | TNFa)
% TNFa(t+1)   = IRF5 & ! IL10s
% TGFb(t+1)   = STAT3 & ! TNFa
% HIF1A(t+1)  = (STAT3 | IL10s) & ! STAT1
% RAGE(t+1)   = HMGB1
%% Initial conditions
x0=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
IFNg = x0(:,1);
GMCSF = x0(:,2);
IL1 = x0(:,3);
LPS = x0(:,4);
IC = x0(:,5);
IL13 = x0(:,6);
IL4 = x0(:,7);
IL10 = x0(:,8);
MCSF = x0(:,9);
HMGB1 = x0(:,10);
IFNgR = x0(:,11);
CSF2Ra = x0(:,12); 
IL1R = x0(:,13);
TLR4 = x0(:,14);
FCgR = x0(:,15);
IL4Ra = x0(:,16);
IL10r = x0(:,17);
MCSFr = x0(:,18);
STAT1 = x0(:,19);
STAT5 = x0(:,20);
NFkb = x0(:,21);
PPARg = x0(:,22);
STAT6 = x0(:,23);
JMJD3 = x0(:,24);
STAT3 = x0(:,25);
IRF3 = x0(:,26);
ERK = x0(:,27);
KLF4 = x0(:,28);
SOCS1 = x0(:,29);
IRF4 = x0(:,30);
IRF5 = x0(:,31);
IL1b = x0(:,32);
IFNab = x0(:,33);
EGF = x0(:,34);
IL12 = x0(:,35);
IL10s = x0(:,36);
TNFa = x0(:,37);
TGFb = x0(:,38);
HIF1A = x0(:,39);
RAGE = x0(:,40);
%% Parameters of simulation
i = 1;                  % Index of simulation
s = 7;                 % number of step
IFNg(:,i) = IFNg;
GMCSF(:,i) = GMCSF;
IL1(:,i) = IL1;
LPS(:,i) = LPS;
IC(:,i) = IC;
IL13(:,i) = IL13;
IL4(:,i) = IL4;
IL10(:,i) = IL10;
MCSF(:,i) = MCSF;
HMGB1(:,i) = HMGB1;
IFNgR(:,i) = IFNgR;
CSF2Ra(:,i) = CSF2Ra; 
IL1R(:,i) = IL1R;
TLR4(:,i) = TLR4;
FCgR(:,i) = FCgR;
IL4Ra(:,i) = IL4Ra;
IL10r(:,i) = IL10r;
MCSFr(:,i) = MCSFr;
STAT1(:,i) = STAT1;
STAT5(:,i) = STAT5;
NFkb(:,i) = NFkb;
PPARg(:,i) = PPARg;
STAT6(:,i) = STAT6;
JMJD3(:,i) = JMJD3;
STAT3(:,i) = STAT3;
IRF3(:,i) = IRF3;
ERK(:,i) = ERK;
KLF4(:,i) = KLF4;
SOCS1(:,i) = SOCS1;
IRF4(:,i) = IRF4;
IRF5(:,i) = IRF5;
IL1b(:,i) = IL1b;
IFNab(:,i) = IFNab;
EGF(:,i) = EGF;
IL12(:,i) = IL12;
IL10s(:,i) = IL10s;
TNFa(:,i) = TNFa;
TGFb(:,i) = TGFb;
HIF1A(:,i) = HIF1A;
RAGE(:,i) = RAGE;
xs = [IFNg(:,1), GMCSF(:,1), IL1(:,1), LPS(:,1), IC(:,1), IL13(:,1),...
    IL4(:,1), IL10(:,1), MCSF(:,1), HMGB1(:,1), IFNgR(:,1), CSF2Ra(:,1),...
    IL1R(:,1), TLR4(:,1), FCgR(:,1), IL4Ra(:,1), IL10r(:,1), MCSFr(:,1),...
    STAT1(:,1), STAT5(:,1), NFkb(:,1), PPARg(:,1), STAT6(:,1), JMJD3(:,1),...
    STAT3(:,1), IRF3(:,1), ERK(:,1), KLF4(:,1), SOCS1(:,1), IRF4(:,1),...
    IRF5(:,1), IL1b(:,1), IFNab(:,1), EGF(:,1), IL12(:,1), IL10s(:,1),...
    TNFa(:,1),TGFb(:,1), HIF1A(:,1), RAGE(:,1)];
%% Boolean simulation of the network
for i = 1:s
    IFNg(:,i+1) = IFNg(:,i);
    GMCSF(:,i+1) = GMCSF(:,i);
    IL1(:,i+1) = IL1(:,i);
    LPS(:,i+1) = LPS(:,i);
    IC(:,i+1) = IC(:,i);
    IL13(:,i+1) = IL13(:,i);
    IL4(:,i+1) = IL4(:,i);
    IL10(:,i+1) = IL10(:,i);
    MCSF(:,i+1) = MCSF(:,i);
    HMGB1(:,i+1) = HMGB1(:,i);
    IFNgR(:,i+1)  = IFNg(:,i) | IFNab(:,i) & not(SOCS1(:,i));
    CSF2Ra(:,i+1) = GMCSF(:,i);
    IL1R(:,i+1)   = IL1(:,i) | IL1b(:,i);
    TLR4(:,i+1)   = LPS(:,i) & not(FCgR(:,i));
    FCgR(:,i+1)   = IC(:,i) & (LPS(:,i) | IL1(:,i));
    IL4Ra(:,i+1)  = IL4(:,i) & IL13(:,i);
    IL10r(:,i+1)  = IL10(:,i) | IL10s(:,i);
    MCSFr(:,i+1)  = MCSF(:,i);
    STAT1(:,i+1)  = IFNgR(:,i) | STAT1(:,i) & not(STAT6(:,i));
    STAT5(:,i+1)  = CSF2Ra(:,i) & not(STAT3(:,i) | IRF4(:,i));
    NFkb(:,i+1)   = (STAT1(:,i) | TNFa(:,i) | TLR4(:,i) | IL1R(:,i)) & not(STAT6(:,i) | FCgR(:,i) | PPARg(:,i) | KLF4(:,i));
    PPARg(:,i+1)  = IL4Ra(:,i) | MCSFr(:,i) | ERK(:,i) & not(STAT6(:,i));
    STAT6(:,i+1)  = IL4Ra(:,i) | MCSFr(:,i);
    JMJD3(:,i+1)  = IL4Ra(:,i) | MCSFr(:,i);
    STAT3(:,i+1)  = (IL10r(:,i) | EGF(:,i) | STAT3(:,i)) & not(FCgR(:,i) | PPARg(:,i));
    IRF3(:,i+1)   = TLR4(:,i);
    ERK(:,i+1)    = FCgR(:,i);
    KLF4(:,i+1)   = STAT6(:,i);
    SOCS1(:,i+1)  = STAT6(:,i) | STAT1(:,i);
    IRF4(:,i+1)   = JMJD3(:,i);
    IRF5(:,i+1)   = STAT5(:,i) & not(IRF4(:,i));
    IL1b(:,i+1)   = NFkb(:,i) | TNFa(:,i);
    IFNab(:,i+1)  = IRF3(:,i);
    EGF(:,i+1)    = ERK(:,i) | STAT3(:,i);
    IL12(:,i+1)   = STAT1(:,i) | STAT5(:,i) | NFkb(:,i);
    IL10s(:,i+1)  = (PPARg(:,i) | STAT3(:,i)) & not(IRF5(:,i) | TNFa(:,i));
    TNFa(:,i+1)   = IRF5(:,i) & not(IL10s(:,i));
    TGFb(:,i+1)   = STAT3(:,i) & not(TNFa(:,i));
    HIF1A(:,i+1)  = (STAT3(:,i) | IL10s(:,i)) & not(STAT1(:,i));
    RAGE(:,i+1)   = HMGB1(:,i);
    xt(i,:)=[IFNg(:,i+1), GMCSF(:,i+1), IL1(:,i+1), LPS(:,i+1), IC(:,i+1),...
        IL13(:,i+1), IL4(:,i+1), IL10(:,i+1), MCSF(:,i+1), HMGB1(:,i+1),...
        IFNgR(:,i+1), CSF2Ra(:,i+1), IL1R(:,i+1), TLR4(:,i+1), FCgR(:,i+1),...
        IL4Ra(:,i+1), IL10r(:,i+1), MCSFr(:,i+1), STAT1(:,i+1), STAT5(:,i+1),...
        NFkb(:,i+1), PPARg(:,i+1), STAT6(:,i+1), JMJD3(:,i+1), STAT3(:,i+1),...
        IRF3(:,i+1), ERK(:,i+1), KLF4(:,i+1), SOCS1(:,i+1), IRF4(:,i+1),...
        IRF5(:,i+1), IL1b(:,i+1), IFNab(:,i+1), EGF(:,i+1), IL12(:,i+1),...
        IL10s(:,i+1), TNFa(:,i+1),TGFb(:,i+1), HIF1A(:,i+1), RAGE(:,i+1)];
end
%% Visualization of dyynamic transition of the network
xt
% Binario a decimal
r = bin2dec(num2str([x0;xt]))
% Decimal a binario 
% dec2bin(6,3)
stairs(r)