%%% ME 599: Cohesive Networks
%%% William Howe
%%% HW 2

close all; clear all; clc;

global n A Ad Bd Zd
%% Problem 1a %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

Ts = max(Set_Time)';

disp('System Deviation Without DSR')
[delta] = cohesiveness(Z,t,Ts,n^2)


%% Problem 1b %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha = 4/Ts;

[K,Bhat]=network_matrix(n,1); 

eK = eig(K);
beta = 0.9*(2*max(1./eK));


tau = Ts/100; 

A = -alpha*beta*K+(1/tau)*(eye(n^2)-beta*K);
Ad = -(1/tau)*(eye(n^2)-beta*K);
Bd = alpha*beta*Bhat;


history = [zeros(n^2,1)];
%time = linspace(0,Tmax,tsteps); 

sol = ddensd(@delay_fun,tau,tau,history,[0,max(t)]);
y = deval(sol,t);
%% Problem 1c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
hold on
plot(t,y')



for i = 1:9
    a = figure(1).Children(2).Children; 
    a(i).Color = [a(i+9).Color];
    a(i).DisplayName = [a(i+9).DisplayName ' With DSR'];
    a(i).LineStyle = '--'; 
    a(i).LineWidth = 2;
    a(i+9).LineWidth = 2;
end
%% Problem 1d %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('System Deviation With DSR')
[delta] = cohesiveness(y',t,Ts,n^2)
%% Bonus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_animation(t,Z,y',0)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% delay function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Tdot] = delay_fun(t,T,Tdelay,Tpdelay)

global n A Ad Bd Zd

% define current and delayed states 
T = T(1:n^2,1);
Tpdelay = Tpdelay(1:n^2,1);

Tdelay = Tdelay(1:n^2,1);

Tdot = A*T+Ad*Tdelay+Bd*Zd;

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
%% Function to animate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data_animation(T,y,z,save_movie)
%%% modified from http://web.mit.edu/8.13/matlab/MatlabTraining_IAP_2012/AGV/DemoFiles/ScriptFiles/html/Part3_Animation.html
f = figure(2);
f.Position = [100 100 1000 800];
ht = sgtitle(sprintf('Time: %0.2f sec', T(1)),'fontsize',30);

subplot(1,2,1)

hh1(1) = line(1, y(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(2) = line(2, y(1,2), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(3) = line(3, y(1,3), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(4) = line(4, y(1,4), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(5) = line(5, y(1,5), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(6) = line(6, y(1,6), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(7) = line(7, y(1,7), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(8) = line(8, y(1,8), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
hh1(9) = line(9, y(1,9), 'Marker', '.', 'MarkerSize', 20, 'Color', 'b');
axis([0 10 0 1]);
title('Without DSR','fontsize',24)
ylabel('Postion','fontsize',24);
xlabel('Element Number','fontsize',24);
xticks([1,2,3,4,5,6,7,8,9]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
b = get(gca,'YTickLabel');
set(gca,'YTickLabel',b,'fontsize',20)


subplot(1,2,2)
hh2(1) = line(1, z(1,1), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
hh2(2) = line(2, z(1,2), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
hh2(3) = line(3, z(1,3), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
hh2(4) = line(4, z(1,4), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
hh2(5) = line(5, z(1,5), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
hh2(6) = line(6, z(1,6), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
hh2(7) = line(7, z(1,7), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
hh2(8) = line(8, z(1,8), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
hh2(9) = line(9, z(1,9), 'Marker', '.', 'MarkerSize', 20, 'Color', [0 .5 0]);
axis([0 10 0 1]);
title('With DSR','fontsize',24)
xlabel('Element Number','fontsize',24);
xticks([1,2,3,4,5,6,7,8,9]);
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',20)
b = get(gca,'YTickLabel');
set(gca,'YTickLabel',b,'fontsize',20)

  
fprintf('Animation (Smart update): %0.2f sec\n', toc);

% Preallocate movie structure.
mov(1:length(T)) = struct('cdata', [], 'colormap', []);

% Animate and add animation frame to the movie structure
for id = 1:length(T)
   % Update XData and YData
   set(hh1(1), 'XData', 1         , 'YData', y(id, 1));
   set(hh1(2), 'XData', 2          , 'YData', y(id, 2));
   set(hh1(3), 'XData', 3          , 'YData', y(id, 3));
   set(hh1(4), 'XData', 4         , 'YData', y(id, 4));
   set(hh1(5), 'XData', 5          , 'YData', y(id, 5));
   set(hh1(6), 'XData', 6          , 'YData', y(id, 6));
   set(hh1(7), 'XData', 7          , 'YData', y(id, 7));
   set(hh1(8), 'XData', 8          , 'YData', y(id, 8));
   set(hh1(9), 'XData', 9          , 'YData', y(id, 9));
   
   set(hh2(1), 'XData', 1         , 'YData', z(id, 1));
   set(hh2(2), 'XData', 2          , 'YData', z(id, 2));
   set(hh2(3), 'XData', 3          , 'YData', z(id, 3));
   set(hh2(4), 'XData', 4         , 'YData', z(id, 4));
   set(hh2(5), 'XData', 5          , 'YData', z(id, 5));
   set(hh2(6), 'XData', 6          , 'YData', z(id, 6));
   set(hh2(7), 'XData', 7          , 'YData', z(id, 7));
   set(hh2(8), 'XData', 8          , 'YData', z(id, 8));
   set(hh2(9), 'XData', 9          , 'YData', z(id, 9));
   set(ht, 'String', sprintf('Time: %0.2f sec', T(id)));
    
   drawnow
   % Get frame as an image
   mov(id) = getframe(gcf);
end

% Create AVI file
if save_movie == 1
    vw = VideoWriter('hw2_bonus_animation.avi');  %taking a guess that you intend to modify the filename each time you write a video
    open(vw);
    writeVideo(vw, mov);
    close(vw);
    clear mov
end
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

