%
% % For Cohesive Networks Class
clear all
% close all
nfig=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%% 
gamma =10;  
N = 10; % number of agents
K = zeros(N,N); % pinned Laplacian
for jj=1:N
    if jj ==1
        K(jj,jj) = 2/2; K(jj,jj+1)=-1/2;
    elseif jj==N
        K(jj,jj) = 1; K(jj,jj-1)=-1;
    else
        K(jj,jj)=2/2; K(jj,jj-1)=-1/2; K(jj,jj+1)=-1/2;
    end
end
Bhat = zeros(N,1); Bhat(1,1) =1/2; 
K = gamma*K;Bhat = Bhat*gamma;
A = -K; 
SysPoles = eig(A)
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % simulate the system with DSR 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% set Sim_start to 0 for using stored data
% set Sim_start to 1 for recomputing
Sim_start =1; 
if Sim_start ==1
    
    
        tmax = 4/abs(max(SysPoles))
        eigK = eig(K)
        delt = tmax/100000; 
        t = 0:delt:2*tmax;  
        % select the control parameters 
        alpha = 4/tmax;             % keep the settling time the same
        tau = delt;           % delay used in simulations 
        beta = 0.9*(4/max(eigK))/(2+alpha*tau);  % selection of beta
        vel = 1;
%return
      % System simulation 
      
% system with delay
Ahat = (eye(size(K))) -delt*alpha*beta*K +(eye(size(K)) -beta*K);
Ahatd = -(eye(size(K)) -beta*K);
Bhatd  = delt*alpha*beta*Bhat;
Astd = (eye(size(K))) - delt*K;
Bstd = delt*Bhat;
ns = length(t);
Is = (3/4)*(2*pi)*(t<=tmax) +(1.25)*(2*pi)*(t>tmax);
Idsr = (pi/2)*ones(N,ns);
Istd = (pi/2)*ones(N,ns);
% initializing the positions 
Xdsr = zeros(N,ns); Ydsr = zeros(N,ns);
Xstd = zeros(N,ns); Ystd = zeros(N,ns);
Xdsr(:,1) = [1:N]'; Xdsr(:,2)=Xdsr(:,1);
Xstd(:,1) = [1:N]'; Xstd(:,2)=Xstd(:,1);
for m=2:(length(t)-1)
    Istd(:,m+1) = Astd*Istd(:,m) + Bstd*Is(m);
    Idsr(:,m+1) = Ahat*Idsr(:,m) + Ahatd*Idsr(:,m-1) + Bhatd*Is(m);
    Xstd(:,m+1) = Xstd(:,m) + delt*vel*cos(Istd(:,m));
    Ystd(:,m+1) = Ystd(:,m) + delt*vel*sin(Istd(:,m));
    Xdsr(:,m+1) = Xdsr(:,m) + delt*vel*cos(Idsr(:,m));
    Ydsr(:,m+1) = Ydsr(:,m) + delt*vel*sin(Idsr(:,m));
end
save L6_Sim_results 
else
    
load L6_Sim_results 
 
end
    nfig=nfig+1; figure(nfig); clf
    plot(t,Is*180/pi,'g',t,Istd*180/pi,'r',t,Idsr*180/pi,'b','LineWidth',3);    
    xlabel('time'), ylabel('orientation')
    %axis([0 max(t) 0 ymax])
    grid
    set(gca,'FontSize',20)
    pause(0.01)
    %saveas(gcf,'../Images/L4_Z_N1_cohesive_response','epsc')
 
 %return   
% quantifying cohesion
[Cstd] = cohesion(t,Istd,Is,N,delt)%  
[Cdsr] = cohesion(t,Idsr,Is,N,delt)%  
%return
%return
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      % plot video
%      
%      filename = 'movie_L6.avi'
%     DeltT = 500; % coarser resampling of time in video
%     [nfig] = storevideo(t,Istd,t,Idsr,N,DeltT,filename,'b',nfig);
% 
%     
%    return
     
    filename = 'movie_L6_XY.avi'
    DeltT = 1000; % coarser resampling of time in video
    [nfig] = storevideoXY(t,Xstd,Ystd,Xdsr,Ydsr,N,DeltT,filename,'b',nfig);
    return
 
    return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function to plot XY positions video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cohesion_loss_normalized] = cohesion(t,I,Is,N,delt)%  
% plot and generate a video   
%
Id = max(max(abs(Is))); % select as max of source value
Ts = max(t); % since we have two steps here
[m,n]= size(I);
Iavg = ones(N,1).*sum(I,1)/N; % average at each time
Idiff = abs(I-Iavg); % difference from average
Cohesion_loss = (sum(sum(Idiff,1))*delt)/Id;
Cohesion_loss_normalized = Cohesion_loss/Ts;        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function to plot XY positions video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nfig] = storevideoXY(t,Xstd,Ystd,Xdsr,Ydsr,N,DeltT,filename,qcolor,nfig)%
% plot and generate a video   
%
           nfig=nfig+1;
           
Xp = max(max(max(Xstd)),max(max(Xdsr)));
Yp = max(max(max(Ystd)),max(max(Ydsr)));
Xn = min(min(min(Xstd)),min(min(Xdsr)));
Yn = min(min(min(Ystd)),min(min(Ydsr)));
%return
           vidObj = VideoWriter(filename);
           open(vidObj);
           for jj=1:DeltT:length(t)
               figure(nfig); clf
               % plot(Xstd(:,jj),Ystd(:,jj),'ro','LineWidth',3);    
                
plot(Xstd(:,jj),Ystd(:,jj),'ro',Xdsr(:,jj),Ydsr(:,jj),'k*','LineWidth',3);    
                axis([Xn Xp Yn Yp])
                grid
                set(gca,'FontSize',20)
                pause(0.1)
                currFrame = getframe;
                writeVideo(vidObj,currFrame);
                pause(0.001) 
           end
                % Close the file.
                close(vidObj);  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  function to plot a video
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nfig] = storevideo(t,I,tdsr,Idsr,N,DeltT,filename,qcolor,nfig)%  
% plot and generate a video   
%
           nfig=nfig+1;
           vidObj = VideoWriter(filename);
           open(vidObj);
           agents = 1:1:N; agents= agents';
           Imax = max(max(I));
           for jj=1:DeltT:length(t)
               figure(nfig); clf
                plot(agents,I(:,jj),'ro',agents,Idsr(:,jj),'ko','LineWidth',3);    
                xlabel('Agents'), ylabel('I')
                axis([0 (max(agents)) 0 Imax])
                grid
                set(gca,'FontSize',20)
                pause(0.1)
                currFrame = getframe;
                writeVideo(vidObj,currFrame);
                pause(0.001) 
           end
                % Close the file.
                close(vidObj);  
end
    
  
      
      
      
      