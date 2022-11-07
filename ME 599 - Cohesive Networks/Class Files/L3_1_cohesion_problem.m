%
% For Cohesive Networks Class
clear all
% close all
nfig=0;

% Figure 1

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


%return
I0 = 0:1:N-1; %I0 = I0*0.5;

Sys = ss(A,B,C,D); 
SysPoles = eig(A)

steady_state = -inv(A)*B

%return
% Part  1 gamma=10 trajectory 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmax = 6/abs(max(SysPoles))
delt = tmax/100; 

t = 0:delt:tmax; 
    
Is = ones(size(t))';
 
[y,t,I]  = lsim(Sys,Is,t);


told = t;
nfig=nfig+1; figure(nfig); clf
%plot(t,I(:,4),'r',t,I(:,3),' r:',t,I(:,2),'g:',t,I(:,1),'g',t,yd,'b','LineWidth',3);
plot(t,Is,'b',t,I(:,1),'g',t,I(:,2),' g:',t,I(:,3),'r:',t,I(:,4),'r','LineWidth',3);
xlabel('time'), ylabel('Position')
%legend('y_4','y_3','y_2','y_1','y_s','location','east')
legend('I_s','I_1','I_2','I_3','I_4','location','east')
%axis([0 max(t) 0 ymax])
grid
title('\gamma =10')
set(gca,'FontSize',20)
pause(0.01)


%return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% faster motion
gamma =100; 
A = -K*gamma; B = Bhat*gamma;


Sys = ss(A,B,C,D); 
delt = 1/100; 
Syspoles = eig(A)


tmax = 4/abs(max(Syspoles))
delt = tmax/100; 
t = 0:delt:tmax;   
yd = ones(size(t))';

 
[y,t,I]  = lsim(Sys,yd,t);

nfig=nfig+1; figure(nfig); clf
plot(t,Is,'b',t,I(:,1),'g',t,I(:,2),' g:',t,I(:,3),'r:',t,I(:,4),'r','LineWidth',3);
xlabel('time'), ylabel('Position')
legend('I_s','I_1','I_2','I_3','I_4','location','east')
title('\gamma =100')
grid
set(gca,'FontSize',20)
pause(0.01)
%saveas(gcf,'../Images/Z_N1_stepresponse','epsc')

return
