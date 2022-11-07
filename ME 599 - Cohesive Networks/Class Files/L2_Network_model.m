
% % For Cohesive Networks Class

clear all
%close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfig=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% define the network 
K = [1 0 0; -1/2 1 -1/2; -1/2 -1/2 1];
Bhat = [1;0;0];
C = eye(3); D = [];

% continuous-time network  

gamma = 1; % alignment gain: settling time = approx 10 s
%gamma = 10; % alignment gain:  settling time = approx 1 s
%gamma = 1000; % stable/ unstable 

 
 
% define the system 
A = -gamma*K; 
B = gamma*Bhat; 
NetSys = ss(A,B,C,D);
% find the step response 
[Y,T,X] = step(NetSys);

% plot the step response
nfig=nfig+1;figure(nfig) 
plot(T,X,'r','LineWidth',2)
xlabel('continuous time'), ylabel('agent states')
grid
set(gca,'FontSize',20)


%return

% Discrete time system
K = [1 0 0; -1/2 1 -1/2; -1/2 -1/2 1];
Bhat = [1;0;0];
C = eye(3); D = [0];E =eye(3);

% different gain values 
dt = 1/100; % sampling time 
P = eye(3) - dt*gamma*K; 
Bd = dt*gamma*Bhat; 

% find the time using dt 
Tdmax = max(T); Td = 0:dt:Tdmax; 

% Initialize the state vector 
Xd =zeros(3,length(Td)); Xd(:,1) = [0;0;0];
U = 1; % input is zero

% find discrete time solution 
for jj=1:1:(length(Td)-1);
    Xd(:,jj+1) = P*Xd(:,jj) + Bd*U;
end

% plot the discrete time solution 
nfig=nfig+1;figure(nfig) 
plot(Td,Xd,'r','LineWidth',2)
xlabel('discrete time'), ylabel('agent states')
grid
%axis([0 max(T) max([min(min(Xd)),-100]) min([max(max(Xd)),100])])
set(gca,'FontSize',20)
