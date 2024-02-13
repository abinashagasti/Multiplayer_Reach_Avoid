%% Performance Metric
% This piece of code helps in measuring the time taken to compute 
% assignment related calculations using brute force method.

sum_time=0;
num_iter=5;

% assign_mat=perms(1:N);
% assign_mat=assign_mat(:,1:M);
assign_mat=nchoosek(1:N,M);
assign_mat=reshape(assign_mat(:,perms(1:M)),[],M);
Nfac=factorial(N);
NMdiff_fac=factorial(N-M);

for avg_timer=1:num_iter

tic;
tstart=tic;

parfor i=1:Nfac/NMdiff_fac
    for j=1:M
        assign_val(i,j)=a(j,assign_mat(i,j));
    end
end

% for j=1:M
%     assign_val(:,j)=a(j,assign_mat(:,j));
% end

sumarr=sum(assign_val,2);
imax=find(sumarr==max(sumarr));
x=zeros(M,N,length(imax));
for i=1:M
    for k=1:length(imax)
        j=assign_mat(imax(k),i);
        x(i,j,k)=1;
    end
end

time_taken=toc(tstart);
sum_time=sum_time+time_taken;

end

avg_time_taken_brute=sum_time/num_iter