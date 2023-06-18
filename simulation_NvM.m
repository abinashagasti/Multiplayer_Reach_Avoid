%% N Pursuer vs M Evader Simulation
% This code shall use the 1P1E situation in addition to a linear
% programming based role assignment for the pursuing and evading teams.

clear
clc
close all

% Initiate player positions

N=4;
M=3;

xp=-15+20*rand(3,N);xp=xp'; % xp,xe = coordinates of pursuers and evaders
xe=-15+20*rand(3,M);xe=xe'; % are in rows
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

% Example 1
% xp=[1,1,1;-1,2,1.5];
% xe=[-1,-0.5,1.5;0,3,1];
% ve=[1;1];vp=[1.5,1.6];

% Example 2
% xp=[-3.9996   -3.1972    8.3785
%    -0.8747   -1.1501   -0.9163];
% xe=[8.9056   -5.6176    7.6481
%    -9.6025   -3.1647    5.3205];
% vp=[1.8428;2.1188];
% ve=[1.4530;1.0102];

% Examples where changing L changes assignment

% Example 1
% xp=[6.4875    9.6533    4.6050;
%    -3.1225    1.6814   -7.8446;
%     8.1262    7.5931    6.3552];
% xe=[-4.7854    1.8871   -9.5497;
%    -1.4948   -3.7456   -6.7703;
%    -6.4247   -1.5423   -8.1154];
% vp=[2.0985;
%     1.9709;
%     2.1959];
% ve=[1.6999;
%     1.6385;
%     1.0336];

% Example 2
% xp=[-1.5760   -0.6957   -2.1588;
%    -6.6190   -7.1848    1.3228;
%    -8.6514    1.2908    0.7815];
% xe=[2.0453   -4.8873   -2.2868;
%     4.0179   -6.1207  -13.7996;
%     2.3350   -2.3762   -7.8985];
% vp=[2.4970;
%     1.7242;
%     2.1525];
% ve=[1.6050;
%     1.3872;
%     1.1422];

% Example 3
% xp=[-6.7681   -2.9472    0.0104;
%    -3.3293   -3.9641   -3.3286;
%    -4.7636  -13.3481   -0.6086];
% xe=[4.9231   -7.9093    4.4252;
%    -8.0710    2.7309   -5.9061;
%    -6.7315  -10.6454  -12.4869];
% vp=[1.7089;
%     2.2261;
%     2.2829];
% ve=[1.6938;
%     1.0098;
%     1.8432];

% Computing static information for assignment

B=zeros(M,N);alpha=zeros(M,N);
a=zeros(M,N);
for i=1:M
    for j=1:N
        alpha(i,j)=ve(i)/vp(j);
        B(i,j)=norm(xe(i,:))^2-alpha(i,j)^2*norm(xp(j,:))^2;
        if B(i,j)>=0 && alpha(i,j)<=1
            a(i,j)=Value(xp(j,:),xe(i,:),alpha(i,j));
        end
    end
end
L=sum(max(a,[],2));
% L=-1;
for i=1:M
    for j=1:N
        if a(i,j)==0
            a(i,j)=-L-1;
        end
    end
end

a

%% Optimal Assignment

% First, we shall find the optimal assignment, because if in the optimal
% assignment the pursuing team wins, then clearly evaders cannot win under
% any condition if the pursuing team plays optimally. 

% f=a';f=f(:);
% b=ones(M+N,1);
% A=zeros(M+N,M*N);
% for i=N+1:M+N
%     A(i,(i-N-1)*N+1:(i-N)*N)=1;
% end
% for j=1:N:M*N
%     A(1:N,j:j+N-1)=eye(N);
% end
% 
% % Optimal assignment
% 
% if N>=M
%     x=linprog(-f,A(1:N,:),b(1:N),A(N+1:N+M,:),b(N+1:N+M),zeros(M*N,1));
%     x=reshape(x,[N,M])'
% end
% if N<M
%     x=linprog(-f,A(N+1:N+M,:),b(N+1:N+M),A(1:N,:),b(1:N),zeros(M*N,1));
%     x=reshape(x,[N,M])';
% end

x=opt_assgn(a,M,N,"primal")

% This was done to ensure that agents in the lower numbner team are
% assigned some agent of the opposite team. But it is not necessary, as if
% some agent is not assigned then even if it were assigned contribution
% would be zero else the solution is not optimal. So nonassignment means
% capture is not possible. 

%% Game of Kind

% Using the optimal assignment we can determine the winning team. If for
% every evader there exists a pursuer that can capture it, then the
% pursuing team wins, else the evading team wins.

win=0; 
% win=0 indicates the pursuer wins, while win=1 indicates the evader wins
check=a.*x;
if min(check(:))<0
    win=1
else
    win=0
end

%% Dual of game 

% Considering the dual of the optimal assignment LP

y=opt_assgn(a,M,N,"dual")
p=zeros(N,1);
for j=1:N
    if sum(x(:,j))==0
        p(j)=y(j);
    else
        i=find(x(:,j)==1);
        p(j)=y(j)+y(N+i)-y(N+M+i);
    end
end

% Check if payoff vector is in core
% 
% S=dec2bin(0:2^N-1)-'0'; % Denotes all possible subsets of N in binary 
% S=[S,zeros(2^N,1)];
% S=characteristic(S,a,M,N);

%% Game of Degree

plot3(0,0,0,'go',LineWidth=3,MarkerSize=8)
hold on
plot3(xp(:,1),xp(:,2),xp(:,3), 'ro',LineWidth=3,MarkerSize=8)
plot3(xe(:,1),xe(:,2),xe(:,3), 'bo',LineWidth=3,MarkerSize=8)
grid on
xlabel('x');ylabel('y');zlabel('z');

if win==0
    disp("The pursuing team wins!")
    T=15;
    z=zeros(T*100+1,9,M);
    V=zeros(M,1);
    for i=1:M
        j=find(x(i,:)==1);
        V(i)=Value(xp(j,:),xe(i,:),alpha(i,j));
        [t,y]=ode45(@agent_dynamics_speed,0:0.01:T,[xe(i,:)';xp(j,:)';ve(i);vp(j);B(i,j)]);
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
