%
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


%return
I0 = 0:1:N-1; %X0 = X0*0.5;



%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simulate the original system without DSR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

Sys = ss(A,B,C,D); 
SysPoles = eig(A)
tmax = 4/abs(max(SysPoles))
delt = tmax/100; 

t = 0:delt:tmax; 
 
Is = ones(size(t))';
 
[y,t,I]  = lsim(Sys,Is,t);

%told = t;
nfig=nfig+1; figure(nfig); clf
plot(t,Is,'b',t,I(:,1),'g',t,I(:,2),' g:',t,I(:,3),'r:',t,I(:,4),'r','LineWidth',3);
plot(t,Is,'b',t,I,'g','LineWidth',3);
xlabel('time'), ylabel('Position')
legend('I_s','I','location','east')
grid
set(gca,'FontSize',20)
pause(0.01)
%saveas(gcf,'../Images/L4_stepresponse_tracking','epsc')

%return


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simulate the system with DSR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% set Sim_start to 0 for using stored data
% set Sim_start to 1 for recomputing
Sim_start =1; 

if Sim_start ==1
    
    
      % select the control parameters 
      eigK =     eig(-A);       % eigenvalues of K
      beta = 0.9*2*max(1./eigK) % selection of beta
      T_s = max(t);             % settling time based on system without DSR
      alpha = 4/T_s             % keep the settling time the same
      tau = alpha/100           % delay used in simulations 
    
      % System simulation 
      
% system with delay
A = -alpha*beta*K + (1/tau)*(eye(size(K)) -beta*K);
Ad = -(1/tau)*(eye(size(K)) -beta*K);
B  = alpha*beta*Bhat;

%
Tmax = max(t);   % maximum simulation time
Id  =  1;        % desired final value
tsteps = length(t);    % number of timesteps
delayT = tau;    % time delay 
save datadelaysoln N delayT Tmax tsteps Id A Ad B  
%

[Idsr,tdsr] = solution_to_delay_diff;
Idsr = Idsr(1:N,:);

save L4_Sim_results 

else
    
load L4_Sim_results 

    
end
    %nfig =nfig+1; figure(nfig); clf;
    %plot(t,x,'r',tdsr,xdsr,'b')
    %xlabel('time')
    %ylabel('X (blue with DSR), (red wihtout DSR)')


    nfig=nfig+1; figure(nfig); clf
    plot(t,Is,'b',t,I,'r',tdsr,Idsr,'k','LineWidth',3);    
    xlabel('time'), ylabel('Position')
    %axis([0 max(t) 0 ymax])
    grid
    set(gca,'FontSize',20)
    pause(0.01)
    %saveas(gcf,'../Images/L4_Z_N1_cohesive_response','epsc')


    nfig=nfig+1; figure(nfig); clf
    plot(t,Is,'b',tdsr,Idsr,'k','LineWidth',3);    
    %legend('I_s','I_1, I_2, I_3, I_4','location','east')
    xlabel('time'), ylabel('Position')
    %axis([0 max(t) 0 ymax])
    grid
    set(gca,'FontSize',20)
    pause(0.01)
    %saveas(gcf,'../Images/Z_N1_cohesive_response_only_2','epsc')
    

    %return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     % plot video

    filename = 'movie_no_DSR.avi'
    DeltT = 1;
    [nfig] = storevideo(t,I,tdsr,Idsr,N,DeltT,filename,'b',nfig);


    return




 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function to simulate delayed system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y, time] =  solution_to_delay_diff;

%save datadelaysoln N delayT Tmax tsteps Is A Ad Bd  
load datadelaysoln

delyp = delayT;
dely = delayT;

history = [zeros(N,1)];
time = linspace(0,Tmax,tsteps); 
% SOL = ddensd( DDEFUN, DELY, DELYP, HISTORY, TSPAN)
sol = ddensd(@ddefun,dely,delyp,history,[0,Tmax]);
y = deval(sol,time);

end % ddex4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function to define the derivative for DDE solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
function dydt = ddefun(t,y,ydel,ypdel)
% Differential equation function

load datadelaysoln
t
% define current and delayed states 
y = y(1:N,1);
ydotdel = ypdel(1:N,1);
ydel = ydel(1:N,1);

dydt = A*y +B*Id +Ad*ydel;
 
end % ddex4de
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function to plot a video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nfig] = storevideo(t,I,tdsr,Idsr,N,DeltT,filename,qcolor,nfig)%  
% plot and generate a video   
%
           nfig=nfig+1;
           vidObj = VideoWriter(filename);
           open(vidObj);
           agents = 1:1:N;
           Imax = max(max(I));
           for jj=1:DeltT:length(t)
               
               figure(nfig); clf
                plot(agents,I(jj,:),'ro',agents,Idsr(:,jj),'ko','LineWidth',3);    
                xlabel('Agents'), ylabel('I')
                axis([0 (max(agents)+1) 0 Imax])
                grid
                set(gca,'FontSize',20)
                pause(0.1)
                currFrame = getframe;
                writeVideo(vidObj,currFrame);
                pause(0.1) 

           end
                % Close the file.
                close(vidObj);  
end
    
  
      
      
      
      
      