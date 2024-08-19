% clear
clc; clear all; close all;

% Simulation time
t0=200;

% Initial conditions
x_0=10; % pray
z_0=1; % predator

% Solving ODEs
[t,y]=ode45(@funodes,[0 t0],[x_0 z_0]);
% Changing variables for presentation
x=y(:,1); z=y(:,2);

% Plotting results
figure(1);
plot(t,x,t,z,'LineWidth',2);
xlabel('Time','FontSize',20);
ylabel('Population number','FontSize',20);
legend({'Pray', 'Predator'},'FontSize',20);
hold on

% ODE function
function dy = funodes(t,y)
% Changing variables for presentation
x=y(1); 
z=y(2);
% Model parameters
k1=1; k2=0.2; k3=0.05; k4=0.1;
% ODEs
dy=zeros(length(y),1);
dy(1)=k1*x-k2*x*z;
dy(2)=k3*x*z-k4*z;
end
