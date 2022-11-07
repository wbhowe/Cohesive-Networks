
% % For Cohesive Networks Class
clear all
% close all
nfig=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
gamma =10;  
N = 4; % number of agents
K = zeros(N,N); % pinned Laplacian
for jj=1:N
    if jj ==1
        K(jj,jj) = 2/2; K(jj,jj+1)=-1/2;
    elseif jj==N
        K(jj,jj) = 1; K(jj,jj-1)=-1;
    else
        K(jj,jj)=2/2; K(jj,jj-1)=-1/2; K(jj,jj+1)=-1/2
    end
end
Bhat = zeros(N,1); Bhat(1,1) =1/2; 
A = -K*gamma; B = Bhat*gamma;
C = zeros(1,N); C(1,N)=1; D=[]; 
% % define the directed network from first class 
% K = [1 0 0; -1 1 0; 0 -1 1];
% Bhat = [1;0;0];
% A = -K*gamma; B = Bhat*gamma;
% C = eye(3); D = [];
%return
I0 = 0:1:N-1; %X0 = X0*0.5;
Sys = ss(A,B,C,D); 
SysPoles = eig(A)
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmax = 4/abs(max(SysPoles))
delt = tmax/100; 
t = 0:delt:tmax; 
 
Is = ones(size(t))';
 
[y,t,I]  = lsim(Sys,Is,t);
%told = t;
nfig=nfig+1; figure(nfig); clf
plot(t,Is,'b',t,I,'r','LineWidth',3);
xlabel('time'), ylabel('Position')
legend('I_s','I','location','east')
grid
set(gca,'FontSize',20)
pause(0.01)
%saveas(gcf,'../Images/L4_stepresponse_tracking','epsc')
%return
% Check stability of DSR 
% using Lambert W function
      eigK =     eig(K);       % eigenvalues of K
      beta = 0.9*2*max(1./eigK) % selection of beta
      T_s = max(t);             % settling time based on system without DSR
      alpha = 4/T_s             % keep the settling time the same
      tau = alpha/100           % delay used in simulations 
      % unstable case 
% check to see if the resulting system is stable -- use Lambert W function
% In the laplace domain
% dx/dt (t)  =  Lambda_i x(t) + Lamda_d,i x(t-tau)
% solutions are 
% s_i = Lambda_i + (1/tau)*W_khat([tau Lamda_d_i e^(-Lambda_i tau]) 
%
lambdaK_i = eigK;
Lambda_i = -alpha*beta*lambdaK_i + (( 1 -beta*lambdaK_i)/tau); 
Lamda_d_i = - (( 1 -beta*lambdaK_i)/tau); 
%
%return
N_lambert = 20; 
nfig=nfig+1;figure(nfig); clf;
for khat = -N_lambert:1:N_lambert
    s_i = Lambda_i + (1/tau)*lambertw(khat,tau*Lamda_d_i.*exp(-Lambda_i*tau))
    if max(real(s_i)) < 0  
        plot(real(s_i),imag(s_i),'go','LineWidth',3); hold on
        plot(real(s_i),-imag(s_i),'go','LineWidth',3); 
    else
        plot(real(s_i),imag(s_i),'ro','LineWidth',3); hold on
        plot(real(s_i),-imag(s_i),'ro','LineWidth',3);
    end
end
hold off; grid
xlabel('real'), ylabel('imaginary')
set(gca,'FontSize',20)
return
 
 