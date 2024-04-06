function [dxi,k] = B2Diff_CLLGRN(t,x,tu,u)
%% Initial parameters
b = 15;
epsilon = 0.5;
dxi = zeros(26,1);
f_x = zeros(26,1);
%Initial condition
IL1=x(1);
IL10=x(2);
HMGB1=x(3);
IFNg=x(4);
GMCSF=x(5);
LPS=x(6);
IC=x(7);
IL13=x(8);
IL4=x(9);
MCSF=x(10);
IFNgR=x(11);
STAT1=x(12);
STAT5=x(13);
NFkb=x(14);
PPARg=x(15);
STAT6=x(16);
JMJD3=x(17);
STAT3=x(18);
IRF4=x(19);
IRF5=x(20);
EGF=x(21);
IL12=x(22);
TNFa=x(23);
TGFb=x(24);
HIF1A=x(25);
RAGE=x(26);
%% Interpolation between tu of controller and t of the solver since u is a time-dependant paramater
u= interp1(tu,u,t);
%% Transformation of Boolean equations into continuos equations
%p1, p2, p3, p4 y p5 NUNCA SE COMENTAN
%M1 signaling
f_x(1) = IL1;
%f_x(1) = IL1+u-IL1*u;  %OR
%f_x(1) = IL1*u;    %AND

%M2 signaling
f_x(2) = IL10;
%f_x(2) = IL10+u-IL10*u; %OR
%f_x(2) = IL10*u;   %AND

%NLC signaling
f_x(3) = HMGB1;
%f_x(3) = HMGB1+u - HMGB1*u; %OR
%f_x(3) = HMGB1*u;  %AND

%M1 signaling
f_x(4) = IFNg;
%f_x(4) = IFNg+u - IFNg*u;  %OR
%f_x(4) = IFNg*u;   %AND

f_x(5) = GMCSF;
%f_x(5) = GMCSF+u - GMCSF*u; %OR
%f_x(5) = GMCSF*u;  %AND

f_x(6) = LPS;
%f_x(6) = LPS+u - LPS*u;    %OR
%f_x(6) = LPS*u;    %AND

%M2 signaling
f_x(7) = IC;
%f_x(7) = IC+u - IC*u;   %OR
%f_x(7) = IC*u;      %AND

f_x(8) = IL13;
%f_x(8) = IL13+u - IL13*u; %OR
%f_x(8) = IL13*u;    %AND

f_x(9) = IL4;
%f_x(9) = IL4+u - IL4*u;    %OR
%f_x(9) = IL4*u;    %AND

%NLC signaling
f_x(10) = MCSF;
%f_x(10) = MCSF+u - MCSF*u;  %OR
%f_x(10) = MCSF*u;     %AND

p1 = IFNg + ((((LPS * (1 -  (IC * (LPS + IL1 - LPS * IL1))))))*...
    (1 - (STAT6 + STAT1 - STAT6 * STAT1))) - IFNg * ((((LPS * (1 -  (IC * ...
    (LPS + IL1 - LPS * IL1)))))) *  (1 - (STAT6 + STAT1 - STAT6 * STAT1)));
f_x(11) = p1;
%f_x(11) = p1+u - p1*u;     %OR
%f_x(11) = p1*u;     %AND

f_x(12) = IFNgR + (STAT1 * (1 - STAT6)) - IFNgR * (STAT1 * (1 - STAT6));
%f_x(12) = (IFNgR + (STAT1 * (1 - STAT6)) - IFNgR * (STAT1 * (1 - STAT6)))+u...
%    - (IFNgR + (STAT1 * (1 - STAT6)) - IFNgR * (STAT1 * (1 - STAT6)))*u;           %OR
%f_x(12) = (IFNgR + (STAT1 * (1 - STAT6)) - IFNgR * (STAT1 * (1 - STAT6)))*u;       %AND

f_x(13) = GMCSF * (1 - (STAT3 + IRF4 - STAT3 *IRF4));
%f_x(13) = (GMCSF * (1 - (STAT3 + IRF4 - STAT3 *IRF4)))+u - (GMCSF * (1 - ...
%    (STAT3 + IRF4 - STAT3 *IRF4)))*u;                                              %OR
%f_x(13) = (GMCSF * (1 - (STAT3 + IRF4 - STAT3 *IRF4)))*u;                          %AND

p2 = ((STAT1 + TNFa - STAT1 * TNFa) + ((LPS * (1 -  (IC * (LPS + ...
    IL1 - LPS * IL1)))) + (IL1 +  (NFkb + TNFa - NFkb * TNFa) - IL1 * ...
    (NFkb + TNFa - NFkb * TNFa)) - (LPS * (1 -  (IC * (LPS + IL1 - LPS * ...
    IL1)))) * (IL1 +  (NFkb + TNFa - NFkb * TNFa) - IL1 * (NFkb + TNFa - ...
    NFkb * TNFa))) - (STAT1 + TNFa - STAT1 * TNFa) * ((LPS * (1 -  (IC * ...
    (LPS + IL1 - LPS * IL1)))) + (IL1 +  (NFkb + TNFa - NFkb * TNFa) - IL1 * ...
    (NFkb + TNFa - NFkb * TNFa)) - (LPS * (1 -  (IC * (LPS + IL1 - LPS * ...
    IL1)))) * (IL1 +  (NFkb + TNFa - NFkb * TNFa) - IL1 * (NFkb + TNFa - NFkb * ...
    TNFa)))) * (1 -  ((STAT6 + (IC * (LPS + IL1 - LPS * IL1)) - STAT6 * (IC * ...
    (LPS + IL1 - LPS * IL1))) + (PPARg + STAT6 - PPARg * STAT6) - (STAT6 + (IC * ...
    (LPS + IL1 - LPS * IL1)) - STAT6 * (IC * (LPS + IL1 - LPS * IL1))) * ...
    (PPARg + STAT6 - PPARg * STAT6)));

f_x(14) = p2;
%f_x(14) = p2+u - p2*u;              %OR                                                                      
%f_x(14) = p2*u;                     %AND

p3 = (IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF + ((IL4 * IL13) + MCSF - ...
    (IL4 * IL13) * MCSF ) + (IC * (LPS + IL1 - LPS * IL1)) - ((IL4 * IL13) + ... 
    MCSF - (IL4 * IL13) * MCSF ) * (IC * (LPS + IL1 - LPS * IL1)) * (1 - STAT6);
f_x(15) = p3;
%f_x(15) = p3+u - p3*u;    %OR
%f_x(15) = p3*u;    %AND

f_x(16) = (IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF;
%f_x(16) = ((IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF)+u - ((IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF)*u;  %OR
%f_x(16) = ((IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF)*u;                                                  %AND

f_x(17) = (IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF;
%f_x(17) = ((IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF)+u - ((IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF)*u;  %OR
%f_x(17) = ((IL4 * IL13) + MCSF - (IL4 * IL13) * MCSF)*u;                                                  %AND

p4 = ((IL10 + ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * ...
    TNFa))) - IL10 * ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * ...
    TNFa)))) + EGF - (IL10 + ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + ...
    TNFa - IRF5 * TNFa))) - IL10 * ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + ...
    TNFa - IRF5 * TNFa)))) * EGF + ((IL10 + ((PPARg + STAT3 - PPARg * STAT3) * ...
    (1 - (IRF5 + TNFa - IRF5 * TNFa))) - IL10 * ((PPARg + STAT3 - PPARg * STAT3) * ...
    (1 - (IRF5 + TNFa - IRF5 * TNFa)))) + EGF - (IL10 + ((PPARg + STAT3 - PPARg * STAT3) * ...
    (1 - (IRF5 + TNFa - IRF5 * TNFa))) - IL10 * ((PPARg + STAT3 - PPARg * STAT3) * ...
    (1 - (IRF5 + TNFa - IRF5 * TNFa)))) * EGF) + STAT3 - ((IL10 + ((PPARg + ...
    STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa))) - IL10 * ((PPARg + ...
    STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa)))) + EGF - (IL10 + ...
    ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa))) - IL10 * ...
    ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa)))) * EGF) * ...
    STAT3 ) * (1 - ((IC * (LPS + IL1 - LPS * IL1)) + PPARg - (IC * (LPS + ...
    IL1 - LPS * IL1)) * PPARg));
f_x(18) = p4;
%f_x(18) = p4+u - p4*u;     %OR
%f_x(18) = p4*u;            %AND

f_x(19) = JMJD3;
%f_x(19) = JMJD3+u - JMJD3*u;   %OR
%f_x(19) = JMJD3*u;             %AND

f_x(20) = STAT5 * (1 - IRF4);
%f_x(20) = (STAT5 * (1 - IRF4))+u - (STAT5 * (1 - IRF4))*u;     %OR
%f_x(20) = (STAT5 * (1 - IRF4))*u;                              %AND

f_x(21) = (IC * (LPS + IL1 - LPS * IL1)) + STAT3 - (IC * (LPS + IL1 - LPS * IL1)) * STAT3;
%f_x(21) = ((IC * (LPS + IL1 - LPS * IL1)) + STAT3 - (IC * (LPS + IL1 - LPS * IL1)) * STAT3)+u...
%    - ((IC * (LPS + IL1 - LPS * IL1)) + STAT3 - (IC * (LPS + IL1 - LPS * IL1)) * STAT3)*u;     %OR
%f_x(21) = ((IC * (LPS + IL1 - LPS * IL1)) + STAT3 - (IC * (LPS + IL1 - LPS * IL1)) * STAT3)*u; %AND

f_x(22) = STAT1 + STAT5 - STAT1 * STAT5 +  (STAT1 + STAT5 - STAT1 * STAT5) + NFkb - (STAT1 + STAT5 - STAT1 * STAT5) * NFkb;
%f_x(22) = (STAT1 + STAT5 - STAT1 * STAT5 +  (STAT1 + STAT5 - STAT1 * STAT5) + NFkb - (STAT1 + STAT5 - STAT1 * STAT5) * NFkb)+u...
%    - (STAT1 + STAT5 - STAT1 * STAT5 +  (STAT1 + STAT5 - STAT1 * STAT5) + NFkb - (STAT1 + STAT5 - STAT1 * STAT5) * NFkb)*u;       %OR
%f_x(22) = (STAT1 + STAT5 - STAT1 * STAT5 +  (STAT1 + STAT5 - STAT1 * STAT5) + NFkb - (STAT1 + STAT5 - STAT1 * STAT5) * NFkb)*u;   %AND

f_x(23) = IRF5 * (1 - ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa))));
%f_x(23) = (IRF5 * (1 - ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa)))))+u...
%    (IRF5 * (1 - ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa)))))*u;         %OR
%f_x(23) = (IRF5 * (1 - ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa)))))*u;   %AND

f_x(24) = STAT3 * (1 - TNFa);
%f_x(24) = (STAT3 * (1 - TNFa))+u - (STAT3 * (1 - TNFa))*u;  %OR
%f_x(24) = (STAT3 * (1 - TNFa))*u;                           %AND

p5 = (STAT3 + ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * ...
    TNFa))) - STAT3 * ((PPARg + STAT3 - PPARg * STAT3) * (1 - (IRF5 + TNFa - IRF5 * TNFa)))) * (1 - STAT1);
f_x(25) = p5;
%f_x(25) = p5+u p5*u;    %OR
%f_x(25) = p5*u;         %AND

%f_x(26) = HMGB1;
%f_x(26) = HMGB1+u - HMGB1*u;    %OR
f_x(26) = HMGB1*u;              %AND

%% ODE´s system
dxi(1)=1/(1+exp(-b*(f_x(1)-epsilon)))-x(1);
dxi(2)=1/(1+exp(-b*(f_x(2)-epsilon)))-x(2);
dxi(3)=1/(1+exp(-b*(f_x(3)-epsilon)))-x(3);
dxi(4)=1/(1+exp(-b*(f_x(4)-epsilon)))-x(4);
dxi(5)=1/(1+exp(-b*(f_x(5)-epsilon)))-x(5);
dxi(6)=1/(1+exp(-b*(f_x(6)-epsilon)))-x(6);
dxi(7)=1/(1+exp(-b*(f_x(7)-epsilon)))-x(7);
dxi(8)=1/(1+exp(-b*(f_x(8)-epsilon)))-x(8);
dxi(9)=1/(1+exp(-b*(f_x(9)-epsilon)))-x(9);
dxi(10)=1/(1+exp(-b*(f_x(10)-epsilon)))-x(10);
dxi(11)=1/(1+exp(-b*(f_x(11)-epsilon)))-x(11);
dxi(12)=1/(1+exp(-b*(f_x(12)-epsilon)))-x(12);
dxi(13)=1/(1+exp(-b*(f_x(13)-epsilon)))-x(13);
dxi(14)=1/(1+exp(-b*(f_x(14)-epsilon)))-x(14);
dxi(15)=1/(1+exp(-b*(f_x(15)-epsilon)))-x(15);
dxi(16)=1/(1+exp(-b*(f_x(16)-epsilon)))-x(16);
dxi(17)=1/(1+exp(-b*(f_x(17)-epsilon)))-x(17);
dxi(18)=1/(1+exp(-b*(f_x(18)-epsilon)))-x(18);
dxi(19)=1/(1+exp(-b*(f_x(19)-epsilon)))-x(19);
dxi(20)=1/(1+exp(-b*(f_x(20)-epsilon)))-x(20);
dxi(21)=1/(1+exp(-b*(f_x(21)-epsilon)))-x(21);
dxi(22)=1/(1+exp(-b*(f_x(22)-epsilon)))-x(22);
dxi(23)=1/(1+exp(-b*(f_x(23)-epsilon)))-x(23);
dxi(24)=1/(1+exp(-b*(f_x(24)-epsilon)))-x(24);
dxi(25)=1/(1+exp(-b*(f_x(25)-epsilon)))-x(25);
dxi(26)=1/(1+exp(-b*(f_x(26)-epsilon)))-x(26);