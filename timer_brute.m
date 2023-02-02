%% Performance Metric
% This piece of code helps in measuring the time taken to compute 
% assignment related calculations using brute force method.

tic;
tstart=tic;
assign_mat=perms(1:N);
assign_mat=assign_mat(:,1:M);
Nfac=factorial(N);

parfor i=1:Nfac
    for j=1:M
        assign_val(i,j)=a(j,assign_mat(i,j));
    end
end
sumarr=sum(assign_val,2);
imax=find(sumarr==max(sumarr));
x=zeros(M,N);
for i=1:M
    j=assign_mat(imax,i);
    x(i,j)=1;
end
time_taken=toc(tstart)