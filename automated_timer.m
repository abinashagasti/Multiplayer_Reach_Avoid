%% Automated code for timer
% Runs a loop to automate computing the average time taken for running all
% the sample cases we consider.

clear
clc
close all

% Test cases

Ncases=[3,5,7,8,9,10,11];
Mcases=[3,4,5,8,6,8,7];
lp_time=zeros(length(Mcases),1);
brute_time=zeros(length(Ncases),1);

for test_case=1:length(Ncases)

    % Initialise positions
    N=Ncases(test_case);
    M=Mcases(test_case);
    xp=-10+20*rand(3,N);xp=xp'; % xp,xe = coordinates of pursuers and 
    xe=-10+20*rand(3,M);xe=xe'; % evaders are in rows
    vp=2+rand(N,1);
    ve=1+rand(M,1);
    
    % Compute 1v1 information
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
    
    % Arrange for LP formulation
    f=a';f=f(:);
    b=ones(M+N,1);
    A=zeros(M+N,M*N);
    for i=1:M+N
        if i>N
            A(i,(i-N-1)*N+1:(i-N)*N)=1;
        end
    end
    for j=1:N:M*N
        A(1:N,j:j+N-1)=eye(N);
    end

    % LP Assignment
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
    
    lp_time(test_case)=sum_time/10;

    % Brute force assignment
    sum_time=0;

    for avg_timer=1:10
    
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
    
    brute_time(test_case)=sum_time/10;

end