%% Performance Metric
% This piece of code helps in measuring the time taken to compute
% coordinates related calculations.

sum_time=0;

for avg_timer=1:10

tstart=tic;

I=zeros(3,M);v=zeros(M);

for i=1:M
    j=find(x(i,:)==1);
    xc=(xe(i,:)-alpha(i,j)^2*xp(j,:))/(1-alpha(i,j)^2);
    rc=(alpha(i,j)/(1-alpha(i,j)^2))*norm(xp(j,:)-xe(i,:));
    I(:,i)=(1-rc/norm(xc))*xc;
    if alpha(i,j)==1
        v(i)=0.5*((norm(xe)^2-norm(xp)^2)/norm(xp-xe));
    else
        v(i)=norm(xc)-rc;
    end
end

time_taken=toc(tstart);
sum_time=sum_time+time_taken;

end

avg_time_taken_coord=sum_time/10