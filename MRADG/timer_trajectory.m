%% Performance Metric
% This piece of code helps in measuring the time taken to compute solution
% trajectories.

win=0; 
check=a.*x;
if min(check(:))<0
    win=1;
else
    win=0;
end

% plot3(0,0,0,'go',LineWidth=3,MarkerSize=8)
% hold on
% plot3(xp(:,1),xp(:,2),xp(:,3), 'ro',LineWidth=3,MarkerSize=8)
% plot3(xe(:,1),xe(:,2),xe(:,3), 'bo',LineWidth=3,MarkerSize=8)
% grid on

if win==0
    disp("The pursuing team wins!")
    z=zeros(1501,8,M);
    for i=1:M
        j=find(x(i,:)==1);
        [t,y]=ode45(@agent_dynamics,[0:0.01:15],[xe(i,:)';xp(j,:)';alpha(i,j);B(i,j)]);
        z(:,:,i)=y;
%         plot3(z(:,1,i),z(:,2,i),z(:,3,i),'b',LineWidth=2)
%         plot3(z(:,4,i),z(:,5,i),z(:,6,i),'r',LineWidth=2)
    end
else
    disp("The evading team wins.")
end