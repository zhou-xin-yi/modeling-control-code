% clear
clc; clear all; close all;

% Simulation time
t0=50;

% Initial conditions
U_0=10^6; % target
I_0=0; % infected
V_0=10; % virus

% Solving ODEs
[t,y]=ode45(@funodes,[0 t0],[U_0 I_0 V_0]);
% Changing variables for presentation
U=y(:,1); I=y(:,2); V=y(:,3);

% Plotting results
figure(1);
plot(t,U,t,I,'LineWidth',2);
xlabel('Time','FontSize',20);
ylabel('Cell number','FontSize',20);
legend({'U','I'},'FontSize',20);
hold on

figure(2);
plot(t,V,'LineWidth',2);
xlabel('Time','FontSize',20);
ylabel('Viral load','FontSize',20);
legend({'V'},'FontSize',20);
hold on

% ODE function
function dy = funodes(t,y)
% Changing variables for presentation
U=y(1); 
I=y(2);
V=y(3);
% Model parameters
U_0=10^6; du=0.01; beta=10^(-7); delta=2; p=1000; c=24; Su=U_0*du;
% ODEs
dy=zeros(length(y),1);
dy(1)=Su-beta*U*V-du*U;
dy(2)=beta*U*V-delta*I;
dy(3)=p*I-c*V;
end
