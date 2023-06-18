function answer = opt_assgn(a,M,N,option)
% Returns the optimal assignment based on matrix a with M evaders and N
% pursuers

f=a';f=f(:);
b=ones(M+N,1);
A=zeros(M+N,M*N);
for i=N+1:M+N
    A(i,(i-N-1)*N+1:(i-N)*N)=1;
end
for j=1:N:M*N
    A(1:N,j:j+N-1)=eye(N);
end
A1=A(1:N,:);A2=A(N+1:N+M,:);
% Optimal assignment

if option=="primal"
    if N>=M
        x=linprog(-f,A(1:N,:),b(1:N),A(N+1:N+M,:),b(N+1:N+M),zeros(M*N,1));
        x=reshape(x,[N,M])';
    end
    if N<M
        x=linprog(-f,A(N+1:N+M,:),b(N+1:N+M),A(1:N,:),b(1:N),zeros(M*N,1));
        x=reshape(x,[N,M])';
    end
    answer=x;
elseif option=="dual"
    y=linprog([b(1:N);b(N+1:N+M);-b(N+1:N+M)],[-A1',-A2',A2'],-f,[],[],zeros(2*M+N));
    answer=y;
else
    error("Wrong option input to opt_assgn function. Input either 'primal' or 'dual'. ")
end


% x_relax=linprog(-f,A,b,[],[],zeros(M*N,1));
% x_relax=reshape(x_relax,[N,M])';
% 
% sum(sum(x-x_relax))

end