%% Performance Metric
% This piece of code helps in measuring the time taken to compute 
% assignment related calculations.

sum_time=0;

for avg_timer=1:10

tic;
tstart=tic;

if N>=M
    x=linprog(-f,A(1:N,:),b(1:N),A(N+1:N+M,:),b(N+1:N+M),zeros(M*N,1));
    x=reshape(x,[N,M])';
end
if N<M
    x=linprog(-f,A(N+1:N+M,:),b(N+1:N+M),A(1:N,:),b(1:N),zeros(M*N,1));
    x=reshape(x,[N,M])';
end

time_taken=toc(tstart);

sum_time=sum_time+time_taken;

end

avg_time_taken_lp=sum_time/10