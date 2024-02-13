%% Storing the permutation matrix sparsely

Nfac=factorial(N);
assign_mat=sparse(Nfac,N);
assign_mat=sparse(perms([1,zeros(1,N-1)]));
assign_mat=assign_mat(:,1:M);

parfor i=1:Nfac
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


%%

m=3;n=10;
C = nchoosek(0:m:m*(n-1),m);
P = perms(1:m);
for ic = 1:size(C,1)
 for ip = 1:size(P,1)
   M = zeros(m,n);
   M(P(ip,:)+C(ic,:)) = 1;
   % Do whatever you want to do with M
 end
end