%%% ME 599: Cohesive Networks
%%% William Howe
%%% HW 4

close all; clear all; clc;
%addpath('/Users/williamhowe/Desktop/UW')
global N Adsr B1dsr B2dsr Zd tau omega
Zd =1;
%% Part 1.a
gamma = 20;
n = 3;
N = n^2; % n being used for pinned laplacian as nxn graph, N is total number of elements


% create the K and B matrix
[K,B] = pinned_laplacian(n,1,gamma);
disp('Pinned Laplacian')
K
disp('Connectivity Matrix')
B

alpha = 1.5;
% creating the state space matrices for the no dsr case
Anodsr = [zeros(N,N),eye(N,N);-alpha^2*K,-2*alpha*K];
Bnodsr = [zeros(N,1);alpha^2*B];
Cnodsr = zeros(1,2*N);Cnodsr(N)=1;
D = [];
sys = ss(Anodsr,Bnodsr,Cnodsr,D);
[~,t,Xnodsr]=step(sys);

figure(1)
plot(t,Xnodsr(:,1:9),'--','linewidth',2)
%% Part 1.b

beta = 10;
tau = .075;



% creating state space dsr matrices
Adsr = [zeros(N,N),eye(N,N);-alpha^2*beta*K,-2*alpha*beta*K]

B1dsr = [zeros(N,1);alpha^2*beta*B]

B2dsr = [zeros(N,N);eye(N,N)-beta*K]



%% Part 1.c
% calculating epsilon lambda
eps = min(abs((beta*eig(K))./(1-beta*eig(K))));

omega = .1*alpha;
%calculating the righ-hand-side of the stability criteria inequality
ubound = eps^(1/2)*alpha*tau;

%plotting range of omega
w = 1:omega/100:1000*omega;
bode_delay = (abs(omega./(1j*w+omega).*(1-exp(-tau*1j*w))));
max_bode_delay = max(bode_delay);

disp('epsilon^(1/2)*alpha*tau')
RH = ubound
disp('sup omega/(s+omega)(1-exp(-tau s))')
LH = max_bode_delay
figure(2)
semilogx(w,bode_delay,w,ubound*ones(size(w)),'r--','linewidth',2)
legend({'$\frac{\omega}{s+\omega}(1-e^{-\tau s})$','$\epsilon^{1/2}*\alpha*\tau$'},'Interpreter','latex','fontsize',16)
ax = gca;
ax.FontSize = 16; 
xlabel('$\omega$','Interpreter','latex','fontsize',30)
%% Part 1.d

history = zeros(4*N,1);
%time = linspace(0,Tmax,tsteps); 

sol = dde23(@delay_fun2,[tau],history,[0, max(t)],[]);
y = deval(sol,t);

figure(1)
hold on
plot(t,y(2*N+1:3*N,:)','linewidth',2)
ax = gca;

ax.FontSize = 16; 
xlabel('Time (s)','fontsize',30)
ylabel('Response','fontsize',30)

%% Part 2.a
B = [1;0];
K = [2 -1;-1 1];

%% Part 2.b
alpha = 2;
N=2;
Anodsr = [zeros(N,N),eye(N,N),zeros(N,N);zeros(N,N),zeros(N,N),eye(N,N);-alpha^3*K,-3*alpha^2*K,-3*alpha*K];
Bnodsr = [zeros(N*2,1);alpha^3*B];
Cnodsr = zeros(1,3*N);Cnodsr(N)=1;
D = [];
sys = ss(Anodsr,Bnodsr,Cnodsr,D);
[~,t,Xnodsr]=step(sys);

figure(3)
plot(t,Xnodsr(:,1:2),'--','linewidth',2)


Adsr = [zeros(N,N),eye(N,N),zeros(N,N);zeros(N,N),zeros(N,N),eye(N,N);-alpha^3*beta*K,-3*alpha^2*beta*K,-3*alpha*beta*K];

B1dsr = [zeros(2*N,1);alpha^3*beta*B];

B2dsr = [zeros(2*N,N);eye(N,N)-beta*K];



eps = min(abs((beta*eig(K))./(1-beta*eig(K))));

omega = .1*alpha;
ubound = eps^(1/2)*alpha*tau;

w = 1:omega/100:1000*omega;
bode_delay = (abs(omega./(1j*w+omega).*(1-exp(-tau*1j*w))));
max_bode_delay = max(bode_delay);

disp('epsilon^(1/2)*alpha*tau')
RH = ubound
disp('sup omega/(s+omega)(1-exp(-tau s))')
LH = max_bode_delay
figure(4)
semilogx(w,bode_delay,w,ubound*ones(size(w)),'r--','linewidth',2)
legend({'$\frac{\omega}{s+\omega}(1-e^{-\tau s})$','$\epsilon^{1/2}*\alpha*\tau$'},'Interpreter','latex','fontsize',16)
ax = gca;
ax.FontSize = 16; 
xlabel('$\omega$','Interpreter','latex','fontsize',30)


history = zeros(6*N,1);
%time = linspace(0,Tmax,tsteps); 

sol = dde23(@delay_fun3,[tau],history,[0, max(t)],[]);
y = deval(sol,t);

figure(3)
hold on
plot(t,y(1:2,:)','linewidth',2)
ax = gca;

ax.FontSize = 16; 
xlabel('Time (s)','fontsize',30)
ylabel('Response','fontsize',30)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% second order delay function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dZdt] = delay_fun2(t,x,xdelay)

global N Adsr B1dsr B2dsr Zd tau omega

% define current and delayed states 
Z = x(2*N+1:3*N,1);
X = x(2*N+1:4*N,1);

Zdelay = xdelay(2*N+1:3*N,1);

Z1delay= xdelay(1:N,1);

Z1 = x(1:N,1);

Z2 = x(N+1:2*N,1);

Z1dot = -omega*Z1+omega*((Z-Zdelay)/tau);

Z2dot = -omega*Z2+omega*((Z1-Z1delay)/tau);

Zdot = Adsr*X+B1dsr*Zd+B2dsr*Z2dot;

dZdt = [Z1dot;Z2dot;Zdot];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% third order delay function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dZdt] = delay_fun3(t,x,xdelay)

global N Adsr B1dsr B2dsr Zd tau omega

% define current and delayed states 
Z = x(1:N,1);
X = x(1:3*N,1);

Zdelay = xdelay(1:N,1);

Z1delay= xdelay(3*N+1:4*N,1);

Z2delay= xdelay(4*N+1:5*N,1);

Z1 = x(3*N+1:4*N,1);

Z2 = x(4*N+1:5*N,1);

Z3 = x(5*N+1:6*N,1);

Z1dot = -omega*Z1+omega*((Z-Zdelay)/tau);

Z2dot = -omega*Z2+omega*((Z1-Z1delay)/tau);

Z3dot = -omega*Z3+omega*((Z2-Z2delay)/tau);

Zdot = Adsr*X+B1dsr*Zd+B2dsr*Z2dot;

dZdt = [Zdot;Z1dot;Z2dot;Z3dot];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to find K and Bhat matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,Bar] = pinned_laplacian(n,source_connections,argin)
    % where n is the nxn size. n^2 is total number of elements. 
    % source_connections is array of integers that represent what the source
    % connectes to. 
    
    Kmat = zeros(n^2,n^2);
    Bar = zeros(n^2,1);
    
    % create a matrix that numerically lists out how the elements are
    % organized. 
    element_map = reshape(1:n^2,n,n);
    
    for i = 1:n^2
        
        % seeing if the current element, i, is connected to the source
        if isempty(find(source_connections==i))==0
            connected = 1;
            Bar(i,1) = 1;
        else
            connected = 0;
        end
        neighbors = 0;
        [x,y] = find(element_map == i);
        
        % sees if element adjacent to element i exists, and if so, what 
        % number it is. Then adds a -1 in that space of the number column
        % of the K matrix for the ith row
        try
            Kmat(i,element_map(x-1,y)) = -1;
            neighbors = neighbors +1;
        end
        try
            Kmat(i,element_map(x+1,y)) = -1;
            neighbors = neighbors +1;
        end
        try
            Kmat(i,element_map(x,y-1)) = -1;
            neighbors = neighbors +1;
        end
        try
            Kmat(i,element_map(x,y+1)) = -1;
            neighbors = neighbors +1;
        end

        Kmat(i,:) = Kmat(i,:);
        Bar(i,1) = Bar(i,1);
    end
    K = Kmat;
    
    try
        K = argin*K;
        Bar = argin*Bar;
    end
    for i = 1:n^2
        K(i,i) = 0-(sum(K(i,:))-Bar(i));
    end
end