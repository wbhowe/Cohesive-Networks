%%% ME 599: Cohesive Networks
%%% William Howe
%%% HW 1

close all; clear all; clc;
addpath('/Users/williamhowe/Desktop/UW')
%% Problem 1a
gamma = 20; 

K = [1,-1/3,0,-1/3,0,0,0,0,0;...
    -1/3,1,-1/3,0,-1/3,0,0,0,0;...
    0,-1/2,1,0,0,-1/2,0,0,0;...
    -1/3,0,0,1,-1/3,0,-1/3,0,0;...
    0,-1/4,0,-1/4,1,-1/4,0,-1/4,0;...
    0,0,-1/3,0,-1/3,1,0,0,-1/3;...
    0,0,0,-1/2,0,0,1,-1/2,0;...
    0,0,0,0,-1/3,0,-1/3,1,-1/3;...
    0,0,0,0,0,-1/2,0,-1/2,1];

A = -gamma*K;
B_hat = [1/3;0;0;0;0;0;0;0;0];
B = gamma*B_hat; 
disp('Eigen values of A Matrix:')
eig_A = eig(A)


%% Problem 1b
C = eye(9);
D = zeros(9,1);

% continuous time state space model
sys = ss(A,B,C,D);

% time steps and input applied
t = 0:.1:10; u = ones(length(t),1);

[y,t] = lsim(sys,u,t);

figure(1)
plot(t,y(:,1),t,y(:,2),t,y(:,3),t,y(:,4),'--',t,y(:,5),'*',t,y(:,6),t,y(:,7),'--',t,y(:,8),'--',t,y(:,9))

legend({'Element 1','Element 2','Element 3','Element 4','Element 5','Element 6','Element 7','Element 8','Element 9'},'location','SouthEast')
xlabel('Time (seconds)')
ylabel('Response')

S = lsiminfo(y,t);

Set_Time = [S.SettlingTime];
disp('System Settling Time (seconds)')
Ts = max(Set_Time)

%% Problem 1d
% time step discretized system
dt = (Ts/1000);

%calculating gamma max based on minimum eigen value 
gamma_max = 2/(dt*min(abs(eig_A)));

P = eye(9) - (gamma_max/10)*dt*K;
Bd = dt*gamma_max/10*B_hat;

% discrete state space model
sys2 = ss(P,Bd,C,D,dt);
t = 0:dt:10; u = ones(length(t),1);
[Y,T] = lsim(sys2,u,t);

figure(2)

plot(T,Y(:,1),T,Y(:,2),T,Y(:,3),T,Y(:,4),'--',T,Y(:,5),'*',T,Y(:,6),T,Y(:,7),'--',T,Y(:,8),'--',T,Y(:,9))
legend({'Element 1','Element 2','Element 3','Element 4','Element 5','Element 6','Element 7','Element 8','Element 9'},'location','SouthEast')
xlabel('Time (seconds)')
ylabel('Response')

%% Problem 1e

K = [1,-1/3,0,-1/3,0,0,0,0,0;...
    -1/4,1,-1/4,0,-1/4,0,0,0,0;...
    0,-1/3,1,0,0,-1/3,0,0,0;...
    -1/4,0,0,1,-1/4,0,-1/4,0,0;...
    0,-1/5,0,-1/5,1,-1/5,0,-1/5,0;...
    0,0,-1/4,0,-1/4,1,0,0,-1/4;...
    0,0,0,-1/3,0,0,1,-1/3,0;...
    0,0,0,0,-1/4,0,-1/4,1,-1/4;...
    0,0,0,0,0,-1/3,0,-1/3,1];


B_hat = [1/3;1/4;1/3;1/4;1/5;1/4;1/3;1/4;1/3];


A = -gamma*K;
B = gamma*B_hat; 
disp('Eigen values of new A Matrix:')
eig_A = eig(A)

sys = ss(A,B,C,D);
t = 0:.01:10; u = ones(length(t),1);

[y,t] = lsim(sys,u,t);

figure(3)
plot(t,y(:,1),t,y(:,2),t,y(:,3),t,y(:,4),'--',t,y(:,5),'*',t,y(:,6),t,y(:,7),'--',t,y(:,8),'--',t,y(:,9))

legend({'Element 1','Element 2','Element 3','Element 4','Element 5','Element 6','Element 7','Element 8','Element 9'},'location','SouthEast')

xlabel('Time (seconds)')
ylabel('Response')


S = lsiminfo(y,t);

Set_Time = [S.SettlingTime];
disp('System Settling Time (seconds)')
Ts = max(Set_Time)


%figure_update

%% Problem 3.1/2: Bonus
% reiterating K for comparison with generalized K function generator
K = [1,-1/3,0,-1/3,0,0,0,0,0;...
    -1/3,1,-1/3,0,-1/3,0,0,0,0;...
    0,-1/2,1,0,0,-1/2,0,0,0;...
    -1/3,0,0,1,-1/3,0,-1/3,0,0;...
    0,-1/4,0,-1/4,1,-1/4,0,-1/4,0;...
    0,0,-1/3,0,-1/3,1,0,0,-1/3;...
    0,0,0,-1/2,0,0,1,-1/2,0;...
    0,0,0,0,-1/3,0,-1/3,1,-1/3;...
    0,0,0,0,0,-1/2,0,-1/2,1];

% calling function that generates K matrix given size n (for nxn system)
[Knew,Bhat_new] = network_matrix(3,1);
disp('K from NxN simulation code')
Knew
disp('K from simulation minus original K matrix')
Knew - K

%% Problem 3.3
[K,Bhat] = network_matrix(10,1);
gamma = 20; 

A = -gamma*K;
B = gamma*Bhat;
C = eye(100);
D = zeros(100,1);
disp('Pole closest to origin:')
eig_A = eig(A);
-min(abs(eig_A))

sys = ss(A,B,C,D);
t = 0:.01:200; u = ones(length(t),1);

[y,t] = lsim(sys,u,t);

figure(4)
plot(t,y)

xlabel('Time (seconds)')
ylabel('Response')


S = lsiminfo(y,t);

Set_Time = [S.SettlingTime];
disp('System Settling Time (seconds)')
Ts = max(Set_Time)

%% Bonus General function
function [K,Bar] = network_matrix(n,source_connections)
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

        Kmat(i,:) = Kmat(i,:)/(neighbors+connected);
        Bar(i,1) = Bar(i,1)/(neighbors+connected);
    end
    K = Kmat + eye(n^2);
end

