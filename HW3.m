%%% ME 599: Cohesive Networks
%%% William Howe
%%% HW 3

close all; clear all; clc;
addpath('/Users/williamhowe/Desktop/UW')

% Problem 1a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 3;
source_connection = 1;
[K,Bhat]=network_matrix(n,source_connection); 

gamma = 20;

A = -gamma*K;
B = gamma*Bhat;
C = eye(n^2);
D = zeros(n^2,1);


sys = ss(A,B,C,D);

Zd = 1;
t = 0:.01:10; u = Zd*ones(length(t),1);

[Z,t] = lsim(sys,u,t);


figure(1)
plot(t,Z,'linewidth',2)

xlabel('Time (seconds)','fontsize',24)
ylabel('Response','fontsize',24)
legend({'Element 1','Element 2','Element 3','Element 4','Element 5','Element 6','Element 7','Element 8','Element 9'},'location','SouthEast','fontsize',16)
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
b = get(gca,'YTickLabel');
set(gca,'YTickLabel',b,'fontsize',20)

S = lsiminfo(Z,t);

Set_Time = [S.SettlingTime];

Lam = eig(A);
disp('Pole Closest to Imaginary axis')
lam = -min(abs(Lam))
disp('Approximate Settling Time')
Ts = 4/min(abs(Lam))

disp('System Settling Time')
Ts = max(Set_Time)

disp('System Deviation Without DSR')
[delta] = cohesiveness(Z,t,Ts,n^2)
%% Problem 1b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = gamma*K;
Bhat = gamma*Bhat;
disp('value of alpha')
alpha = 4/Ts
disp('sampling time')
dt = Ts/1000
disp('maximum eigen value of K')
lam_k = max(eig(K))
disp('maximum value of beta')
beta_star = 4/(lam_k*(alpha*dt+2))
disp('beta value used for system')
beta = .9*beta_star


Ahat = (eye(n^2)-alpha*beta*K*dt+eye(n^2)-beta*K);
Adhat = -(eye(n^2)-beta*K);
Bdhat = alpha*beta*dt*Bhat;

t = 0:dt:Ts;
Z_dt = zeros(9,length(t));

for i = 2:length(Z_dt)-1
    Z_dt(:,i+1) = Ahat*Z_dt(:,i)+Adhat*Z_dt(:,i-1)+Bdhat*Zd;
    
end
disp('System Deviation With DSR')
[delta] = cohesiveness(Z_dt',t,Ts,n^2)

hold on
plot(t,Z_dt,'--')

for i = 1:9
    a = figure(1).Children(2).Children; 
    a(i).Color = [a(i+9).Color];
    a(i).DisplayName = [a(i+9).DisplayName ' With DSR'];
    a(i).LineStyle = '--'; 
    a(i).LineWidth = 2;
    a(i+9).LineWidth = 2;
end
%% Problem 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation formation of nodes tracking target

n = 1000; % number of steps
s = 50; % offset leader from target
x = linspace(0,10,n); % target x trajectory
y = x.^2; % target y trajecotry

% position target along trajectory
target =[x(s:n);y(s:n)];

% position leader node along trajectory
leader = [x(1:n-s+1);y(1:n-s+1)];


% distance between nodes
d = 2;
% angle of V formation
phi = 60*pi/180;
% number of 
N=4;

% creating initial line objects
figure(2)
axis([-5 20 -5 100])
xlabel('X','fontsize',24)
ylabel('Y','fontsize',24)
title('Simulation of Virtual Node formation','fontsize',30)
hh1 = line(leader(1,1), leader(2,1), 'Marker', 'o', 'MarkerSize', 20, 'Color', 'g');
hh2 = line(target(1,1), target(2,1), 'Marker', '*', 'MarkerSize', 20, 'Color', 'r');
for i = 1:N
    hh3(i) = line(0, 0);
    hh4(i) = line(0, 0);
end

for i = 1:length(target(2,:))
    
    set(hh1, 'XData', leader(1,i), 'YData', leader(2,i));
    set(hh2, 'XData', target(1,i), 'YData', target(2,i));


    [qju,qjn]=virtual_nodes(target(:,i),leader(:,i),phi,d,N);

    for j = 1:N
        set(hh3(j), 'XData', qju(1,j), 'YData', qju(2,j), 'Marker', 'o', 'MarkerSize', 20, 'Color', 'b');
        set(hh4(j), 'XData', qjn(1,j), 'YData', qjn(2,j), 'Marker', 'o', 'MarkerSize', 20, 'Color', 'b');
    end
    
    drawnow
    mov(i) = getframe(gcf);
    
end

save_movie = 0;
if save_movie == 1
    vw = VideoWriter('hw3_pr2_virtual_node_simulation.avi');  %taking a guess that you intend to modify the filename each time you write a video
    open(vw);
    writeVideo(vw, mov);
    close(vw);
    clear mov
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation of element moving to virtual node
x = 1; %linspace(0,10,n);
y = x.^2;
% position targe
target =[x;y];

% position leader node
leader = [x-2;y-2];
d = 2;
phi = 60*pi/180;
N=4;

[qju,qjn]=virtual_nodes(target,leader,phi,d,N);

nodes = [leader,qju,qjn];
element =randi([-5,5],1,1)*5*rand(2,1);
figure(3)

plot(nodes(1,:),nodes(2,:),'o','MarkerSize',10)
xlabel('X','fontsize',24)
ylabel('Y','fontsize',24)
title('Simulation of Virtual Node formation','fontsize',30)
hh1 = line(element(1), element(2), 'Marker', '*', 'MarkerSize', 10, 'Color', 'g');

axis([-25 25 -25 25])
[ind] = closest_node(element,nodes);
x2 = 0;
[el] = element_model(element,nodes(:,ind),x2,.01,1,1);
for i = 1:length(el(1,:))
        set(hh1, 'XData', el(1,i), 'YData', el(2,i));
        drawnow
        mov(i) = getframe(gcf);
end

if save_movie == 1
    vw = VideoWriter('hw3_pr2_element_node_tracking.avi');  %taking a guess that you intend to modify the filename each time you write a video
    open(vw);
    writeVideo(vw, mov);
    close(vw);
    clear mov
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bonus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
save_movie = 0;

x = zeros(length(0:.01:10),9);
y = zeros(length(0:.01:10),9);
x_dsr = zeros(length(0:.01:10),9);
y_dsr = zeros(length(0:.01:10),9);


x(1,:)=[0:8];
x_dsr(1,:)=[0:8];

V = 1;
dt = .01;
T = 0:dt:10;

Z_dsr = Z_dt';

for i = 2: length(T)
    dx = cos(Z(i-1,:))*V;
    dy = sin(Z(i-1,:))*V;
    
    dx_dsr = cos(Z_dsr(i-1,:))*V;
    dy_dsr = sin(Z_dsr(i-1,:))*V;
    
    x(i,:) = dt*dx+x(i-1,:);
    y(i,:) = dt*dy+y(i-1,:);
    
    x_dsr(i,:) = dt*dx_dsr+x_dsr(i-1,:);
    y_dsr(i,:) = dt*dy_dsr+y_dsr(i-1,:);
    
    
end

f=figure(4);
f.Position = [100 100 1000 800];

ht = title(sprintf('Time: %0.2f sec', T(1)),'fontsize',30);
xlabel('X','fontsize',30)
ylabel('Y','fontsize',30)
axis([0,20,0,10])
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
b = get(gca,'YTickLabel');
set(gca,'YTickLabel',b,'fontsize',20)
for i = 1:9
    hh1(i) = line(x(1,i), y(1,i), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
    hh2(i) = line(x_dsr(1,i), y_dsr(1,i), 'Marker', '.', 'MarkerSize', 20, 'Color', 'r');
    hh3(i) = line(x(1,i), y(1,i), 'LineStyle', '--', 'LineWidth', 1, 'Color', 'b');
    hh4(i) = line(x_dsr(1,i), y_dsr(1,i), 'LineStyle','--', 'LineWidth', 1, 'Color', 'r');
end
mov(1:length(length(T))) = struct('cdata', [], 'colormap', []);

for i = 1: length(T)
    for j = 1:9
        
        set(hh1(j), 'XData', x(i,j), 'YData', y(i,j));
        set(hh2(j), 'XData', x_dsr(i,j), 'YData', y_dsr(i,j));
        if i<=100
            g=1:i;
        else
            g=i-100:i;
        end
        set(hh3(j), 'XData', x(g,j), 'YData', y(g,j));
        set(hh4(j), 'XData', x_dsr(g,j), 'YData', y_dsr(g,j));
        if i==1 && j==1
            legend({'Without DSR','With DSR'})
        end
    end
    
    ht = title(sprintf('Time: %0.2f sec', T(i)),'fontsize',30);
    set(ht, 'String', sprintf('Time: %0.2f sec', T(i)));

    drawnow
    % Get frame as an image
    mov(i) = getframe(gcf);
end

% Create AVI file
if save_movie == 1
    vw = VideoWriter('hw3_bonus_animation.avi');  %taking a guess that you intend to modify the filename each time you write a video
    open(vw);
    writeVideo(vw, mov);
    close(vw);
    clear mov
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to determine deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delta] = cohesiveness(T,t,Ts,m)
Tbar = (1/m)*(sum(T'))';

f = sum(abs(T-Tbar*ones(1,m))')';

delta = (1/Ts)*trapz(t,f);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to find K and Bhat matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [K,Bar] = network_matrix(n,source_connections,argin)
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
    
    try
        K = argin*K;
        Bar = argin*Bar;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to move element to closest node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x1] = element_model(x1,node,x2,dt,kp,kv)
    i=2;
    while norm((x1(:,i-1)-node),2)>.01

        x2 = dt*(-kp*(x1(:,i-1)-node)-kv*(x2))+x2;
        x1(:,i) = dt*x2+x1(:,i-1);
        i=i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to determine closest node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ind] = closest_node(element,nodes)

    for i = 1:length(nodes(1,:))
        dnorm(:,i) = norm((element-nodes(:,i)),2);
    end
    [~,ind] = min(dnorm);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function to create virtual nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [qju,qjn]=virtual_nodes(Pt,Pl,phi,d,N)
    theta = atan((Pt(2)-Pl(2))/(Pt(1)-Pl(1)));
    delta = phi/2;

    qu_prime = d*[cos(delta);sin(delta)];
    qn_prime = d*[-cos(delta);sin(delta)];


    R2 = [cos(theta),sin(theta);-sin(theta),cos(theta)];
    R1 = [cos(theta),-sin(theta);sin(theta),cos(theta)];

    qu = Pl+R1*qu_prime;
    qn = Pl+R1*qn_prime;

for j = 1:N
    qju(:,j) = (1+j)*Pl-j*qu;
    qjn(:,j) = (1-j)*Pl+j*qn;
end

end
