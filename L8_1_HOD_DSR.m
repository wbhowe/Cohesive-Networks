% Simulations for problem in Automatica Paper
% Second order Complex Case system

    clear all 
    nfig=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  Step 1: stimulate step response of the system without DSR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    N=6;  nleader=1;
    K = zeros(N,N);     % initialize the K matrix
    B = zeros(N,1);     % initialize the B matrix
       % K(1,1:6)    = [0 0 0 0 0 0];
        K(1,1:6)    = [3 0 0 -1 -1 -1];
        K(2,1:3)  = [-1  2 -1];
        K(3,1:3)  = [-1 -1  2];
        K(4,2:4)  = [-1 0  1]; 
        K(5,3:5)  = [-1 0  1];
        K(6,2:6)  = [-1  -1  -1  -1 4];
    % Pinning the leader 
    K(nleader,nleader) = K(nleader,nleader)+1;    
    B(nleader) = 1;
    
    eig(K)
    %return
    % Need to choose beta to satisfy: \beta  > \frac{1}{\rho_{i} \cos{\phi_i} }, 
    % say 2 times the maximum of rhs over all i 
    lambda = eig(K);                    % the eigenvalues of K
    Ts = 7.5;
    tau   = Ts/100; 
    
    % similar settling time but for the no-DSR second order system
    alpha_no_dsr = 1.69; % selected numerically 
    beta =20;  
    alpha =0.5344;  % alpha_no_dsr/sqrt(beta) 

    %return
    alpha_hat1 =  2*alpha; alpha_hat0 = alpha*alpha;
    alpha_hat1_nodsr =  2*alpha_no_dsr; 
    alpha_hat0_nodsr = alpha_no_dsr*alpha_no_dsr;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%  Step 2: simulate the system without the delayed control scheme Second Order System 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
    % simulate the system without the delayed control scheme
    zK = zeros(size(K));  iK = eye(size(K)); % zero and identity matrices
    % define the system

    Anodsr = [zK iK; -alpha_hat0_nodsr*K -alpha_hat1_nodsr*K]; 
    Bnodsr = [0*B; alpha_hat0_nodsr*B]; Cnodsr = zeros(1,2*N); Cnodsr(1,N)=1;  Dnodsr = []; 
    Cvel = zeros(1,2*N); Cvel(1,N)=1; 
    Sysnodelay = ss(Anodsr, Bnodsr, Cnodsr, Dnodsr);
    eigA = eig(Anodsr);

    % Find the setting time and redo step with specified time vector
    [Y,t,X] = step(Sysnodelay); 
    tsteps = 2000;  % time steps in the simulations 
    % redo step response with specified time vector
    tmax = max(t); t=1:1:tsteps; t = t*tmax/tsteps; 
    [Y,t,X] = step(Sysnodelay,t);
    Z_no_dsr = X(:,1:N);
    Zdot_no_dsr = X(:,N+1:2*N);
    U_no_dsr = -alpha_hat0_nodsr*K*Z_no_dsr' -alpha_hat1_nodsr*K*Zdot_no_dsr' +alpha_hat0_nodsr*B;

    %return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Step 3: Stability of system with DSR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    eigK = eig(K);
    eps_lamdba = min(abs(beta*eigK./(1-beta*eigK)))
    %%% check for stability         
    Omega = .1*alpha;  % filter meets stability condition
    %Omega = 11.3*alpha;  % filter does not meet stability condition
    w  =1:Omega/10:100*Omega; 
    fs = Omega./abs((j*w) + Omega) ;
    lhs = (1-exp(-tau*j*w)).*fs;
    lhs = max(abs(lhs));
    rhs = tau*(eps_lamdba.^(1/2))*alpha;
    if rhs >= lhs
        lhs_rhs = [ lhs rhs]
        disp('stability condition is met')
    else
        lhs_rhs = [ lhs rhs]
        disp('Problem: stability condition is not met')
        %return
    end
    
%return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Step 4: simulate the system with DSR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Adsr  =  [zK iK; -beta*alpha_hat0*K -beta*alpha_hat1*K]; 
    Bdsr  =  [0*B; beta*alpha_hat0*B];  
    B2dsr =  [zK; iK-beta*K]; 
  
    Nsim = N*4;    
    history = [zeros(Nsim,1)];
    Param_dde = [alpha tau beta N Omega];
    %sol = dde23(ddefun,lags,history,tspan,options)
    sol = dde23(@ddefun,[tau],history,[0, tmax],[],Adsr,Bdsr,B2dsr,Param_dde);
    y_dsr = deval(sol,t);
    % Definitions in dde23 function dydt = ddefun(t,y,ydel,Adsr,Bdsr,B2dsr,K,Param)
    %alpha = Param(1); tau = Param(2); beta = Param(3); 
    %N = Param(4); Omega = Param(5); 
    
        
    Z_dsr = y_dsr(1:N,:);  
    Zdot_dsr = y_dsr(N+1:2*N,:);
    Zhat2_dsr = y_dsr(3*N+1:4*N,:);
    U_dsr = -beta*alpha_hat0*K*Z_dsr -beta*alpha_hat1*K*Zdot_dsr +beta*alpha_hat0*B  + (eye(size(K))-beta*K)*Zhat2_dsr;

    
    nfig=nfig+1; figure(nfig)
    plot(t,Z_no_dsr,'r:',t,Z_dsr,'b','LineWidth',4);  
    xlabel('time'); ylabel('Z'); 
    set(gca,'FontSize',20)
    pause(0.01)
    % Blue(DSR) Red(no-DSR)
    %legend('no DSR','with DSR')

    nfig=nfig+1; figure(nfig)
    semilogx(t,U_no_dsr,'r:',t,U_dsr,'b', 'LineWidth',4); 
    xlabel('time'); ylabel('U');
    set(gca,'FontSize',20)
    pause(0.01)
    %Blue(DSR) Red(no-DSR)
    %saveas(gcf,'fig_hdf_U_comparison_rev_resp','epsc')

    Max_input_U_dsr = max(max(U_dsr))
    Max_input_U_no_dsr = max(max(U_no_dsr))
    
    return
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function to define the derivative for DDE solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dZdt = ddefun(t,y,ydel,Adsr,Bdsr,B2dsr,Param)
% Differential equation function
alpha = Param(1); tau = Param(2); beta = Param(3); 
N = Param(4); Omega = Param(5);  

% define current and delayed states X=[P;Pdot], Phat1, Phat2
zs =1;              % source is unit step
X = y(1:N*2,1);     % this is X =[Z(t);Zdot(t)]

Z = y(1:N,1);       % this is Z(t)
Zdel = ydel(1:N,1); % this is Z(t-tau)
Zhat1 = y((N*2 +1):N*3,1);       % this is Zhat1(t)
Zhat1del = ydel((N*2 +1):N*3,1); % this is Zhat1(t-tau)
Zhat2 = y((N*3 +1):N*4,1);       % this is Zhat2(t)
Zhat2del = ydel((N*3 +1):N*4,1); % this is Zhat2(t-tau)

% equations for the two filters
dot_Zhat1 = -Omega*Zhat1 +Omega*(Z -Zdel)/tau; 
dot_Zhat2 = -Omega*Zhat2 +Omega*(Zhat1 -Zhat1del)/tau; 
% equations for the network dynamics d/dt X = d/dt [Z(t);Zdot(t)]
dot_X = Adsr*X +Bdsr*zs +B2dsr*Zhat2;

dZdt = [dot_X ; dot_Zhat1; dot_Zhat2];
end % ddex4de


   