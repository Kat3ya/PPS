clc
clear all
%% Solving ODES TAM - CLL
%i_c = [0 0 0 0 1 0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 0 1 1 0 0 0]; %Initial condition M1
%i_c = [0 0 0 0 0 0 0 1 1 0 0 0 0 0 1 1 1 0 1 0 0 0 0 0 1 0]; %Initial condition M2
%i_c = [0 0 0 0 1 0 0 0 0 0 0 0 0 1 0 0 0 1 0 0 1 1 0 1 1 0]; %Initial condition NLC
%i_c = [1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 1 0 0 0 0]; %Initial condition M0
options = odeset('RelTol',1e-3,'AbsTol',[1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4...
    1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4 1e-4], 'MaxStep',0.01); %Configuration of ODE45 solver
[T,X] = ode45(@(t,x) TAMCLL_simulation(t,x),[0 12],i_c,options); %ODE45 solver 

figure
plot(T,X(:,1))
hold on
plot(T,X(:,2))
plot(T,X(:,3))
plot(T,X(:,4))
plot(T,X(:,5))
plot(T,X(:,6))
plot(T,X(:,7))
plot(T,X(:,8))
plot(T,X(:,9))
plot(T,X(:,10))
plot(T,X(:,11))
plot(T,X(:,12))
plot(T,X(:,13))
plot(T,X(:,14))
plot(T,X(:,15))
plot(T,X(:,16))
plot(T,X(:,17))
plot(T,X(:,18))
plot(T,X(:,19))
plot(T,X(:,20))
plot(T,X(:,21))
plot(T,X(:,22))
plot(T,X(:,23))
plot(T,X(:,24))
plot(T,X(:,25))
plot(T,X(:,26))



    
    
    
    