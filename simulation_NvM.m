%% N Pursuer vs M Evader Simulation
% This code shall use the 1P1E situation in addition to a linear
% programming based role assignment for the pursuing and evading teams.

clear
clc
close all

% Initiate player positions

N=2;
M=2;

xp=-10+20*rand(3,N);xp=xp'; % xp,xe = coordinates of pursuers and evaders
xe=-10+20*rand(3,M);xe=xe'; % are in rows
vp=1.5+rand(N,1);
ve=1+rand(M,1);

% For numerical illustrations

% First standard 3v3 example
% xp=[-8.0258,-4.7626,-3.2929;3.5946,-7.2689,4.4245;-7.8648,3.0751,-0.1165];
% xe=[5.5810,4.3007,8.0744;7.8185,-3.3167,3.9749;-6.0438,-9.3892,4.8815];
% vp=[2;1.9799;2.4047];
% ve=[1.6099;1.6177;1.8594];


% Dispersal surface example
% xp=[1 0 0;1 0 0.5;1 0 -0.5];
% xe=[0.75 1 0;0.75 -1 0];
% vp=ones(3);ve=0.5*ones(2);

% Subset of evaders superior to pursuers with N=5, M=3 example
% xp=[-7.8391   -0.8025   -0.9823;
%     1.0228    6.1081    4.0170;
%     7.4447   -8.9562   -5.6064;
%    -0.8072    9.1707    5.8009;
%    -0.9625   -3.3314   -8.8181   ];
% xe=[4.8181    0.1359   -6.0015;
%    -1.4561   -6.6262    5.0339;
%    -2.6330    8.8364   -9.6565];
% vp=[1.8829;
%     1.5846;
%     2.2339;
%     1.8320;
%     2.3397];
% ve=[1.3717;
%     1.8282;
%     1.1765];

% Subset of evaders superior to pursuers with N=3, M=3 example
% xp=[-7.1447    0.0910    2.2138
%     4.0759   -2.3332    4.5737
%     7.7457   -8.8830   -7.2357];
% xe=[7.2613   -1.5651   -1.7737
%     9.1828    5.0049    9.6199
%    -5.3297   -8.0755   -2.3084];
% vp=[%2.0003
%     1.8567
%     2.0703
%     2.4766];
% ve=[1.4929
%     1.4009
%     1.9950];

% Subset of evaders superior to pursuers with N=3,M=3 second example
% xp=[5.3260    5.0261   -7.2227;
%    -3.0136   -6.9732   -0.0656;
%     6.1730    2.6574    3.7680];
% xe=[2.7914    4.5864    7.1969;
%     2.5391   -6.3882    1.4661;
%    -6.7287    8.1210   -8.4531];
% vp=[%1.8385;
%     1.7546
%     2.0806;
%     1.9752];
% ve=[1.8053;
%     1.5308;
%     1.2273];

% Nonoptimal play
% xp=[1,1,1;-1,2,1.5];
% xe=[-1,-0.5,1.5;0,3,1];
% ve=[1;1];vp=[1.5,1.6];

% Nonoptimal play second example
xp=[-3.9996   -3.1972    8.3785
   -0.8747   -1.1501   -0.9163];
xe=[8.9056   -5.6176    7.6481
   -9.6025   -3.1647    5.3205];
vp=[1.8428;2.1188];
ve=[1.4530;1.0102];

% Compute static information

% Find alpha and mu matrix that denotes if evader i can be 
% captured by pursuer j?
% B=zeros(M,N);mu=zeros(M,N);alpha=zeros(M,N);
% h=zeros(M,N);c=zeros(M,1);
% for i=1:M
%     for j=1:N
%         alpha(i,j)=ve(i)/vp(j);
%         B(i,j)=norm(xe(i,:))^2-alpha(i,j)^2*norm(xp(j,:))^2;
%         if B(i,j)>=0 && alpha(i,j)<=1
%             mu(i,j)=1;
%         else
%             mu(i,j)=0;
%         end
%         if mu(i,j)==1
%             h(i,j)=Value(xp(j,:),xe(i,:),alpha(i,j));
%         else
%             h(i,j)=-100;
%         end
%         end
%         c(i)=min(h(i,:));
% end
% Before being able to find the solution of the game of kind, one has to
% look across all possible assignments of evaders to pusuers, and see if
% there exists any assignment where all evaders can be captured. Simply
% taking the least value of the mu matrix does not tell us whether all
% evaders can be captured or not. For a very simple case, take a situation
% where M>N and all mu values are positive, then clearly their min will
% also be positive, but that does not imply that the game shall be won by
% the pursuers, as clearly at least M-N evaders cannot be captured.

% Alternate for previous section

B=zeros(M,N);alpha=zeros(M,N);
a=zeros(M,N);
for i=1:M
    for j=1:N
        alpha(i,j)=ve(i)/vp(j);
        B(i,j)=norm(xe(i,:))^2-alpha(i,j)^2*norm(xp(j,:))^2;
        if B(i,j)>=0 && alpha(i,j)<=1
            a(i,j)=Value(xp(j,:),xe(i,:),alpha(i,j));
        else
            a(i,j)=-100;
        end
    end
end
a
% Optimal Assignment

% First, we shall find the optimal assignment, because if in the optimal
% assignment the pursuing team wins, then clearly evaders cannot win under
% any condition if the pursuing team plays optimally. 
% for j=1:N
%     a(:,j)=h(:,j)-c;
% end
% a=h;
f=a';f=f(:);
b=ones(M+N,1);
A=zeros(M+N,M*N);
for i=N+1:M+N
    A(i,(i-N-1)*N+1:(i-N)*N)=1;
end
for j=1:N:M*N
    A(1:N,j:j+N-1)=eye(N);
end

% Optimal assignment

if N>=M
    x=linprog(-f,A(1:N,:),b(1:N),A(N+1:N+M,:),b(N+1:N+M),zeros(M*N,1));
    x=reshape(x,[N,M])'
end
if N<M
    x=linprog(-f,A(N+1:N+M,:),b(N+1:N+M),A(1:N,:),b(1:N),zeros(M*N,1));
    x=reshape(x,[N,M])';
end

% This was done to ensure that agents in the lower numbner team are
% assigned some agent of the opposite team. But it is not necessary, as if
% some agent is not assigned then even if it were assigned contribution
% would be zero else the solution is not optimal. So nonassignment means
% capture is not possible. 

x_relax=linprog(-f,A,b,[],[],zeros(M*N,1));
x_relax=reshape(x_relax,[N,M])';
% 
% sum(sum(x-x_relax))

% Game of Kind

% Using the optimal assignment we can determine the winning team. If for
% every evader there exists a pursuer that can capture it, then the
% pursuing team wins, else the evading team wins.

win=0; 
% win=0 indicates the pursuer wins, while win=1 indicates the evader wins
check=a.*x;
if min(check(:))<0
    win=1;
else
    win=0;
end

%% Game of Degree

plot3(0,0,0,'go',LineWidth=3,MarkerSize=8)
hold on
plot3(xp(:,1),xp(:,2),xp(:,3), 'ro',LineWidth=3,MarkerSize=8)
plot3(xe(:,1),xe(:,2),xe(:,3), 'bo',LineWidth=3,MarkerSize=8)
grid on
xlabel('x');ylabel('y');zlabel('z');

if win==0
    disp("The pursuing team wins!")
    z=zeros(1501,9,M);
    V=zeros(M,1);
    for i=1:M
        j=find(x(i,:)==1);
        V(i)=Value(xp(j,:),xe(i,:),alpha(i,j));
        [t,y]=ode45(@agent_dynamics_speed,[0:0.01:15],[xe(i,:)';xp(j,:)';ve(i);vp(j);B(i,j)]);
        z(:,:,i)=y;
        plot3(z(:,1,i),z(:,2,i),z(:,3,i),'b',LineWidth=2)
        plot3(z(:,4,i),z(:,5,i),z(:,6,i),'r',LineWidth=2)
    end
    netValue=sum(V)
else
    disp("The evading team wins.")
end

%% Writing data

% writematrix(z(:,1:3,1),'xe1.csv','Delimiter','tab')
% writematrix(z(:,1:3,2),'xe2.csv','Delimiter','tab')
% writematrix(z(:,1:3,3),'xe3.csv','Delimiter','tab')
% writematrix(z(:,4:6,1),'xp1.csv','Delimiter','tab')
% writematrix(z(:,4:6,2),'xp2.csv','Delimiter','tab')
% writematrix(z(:,4:6,3),'xp3.csv','Delimiter','tab')
