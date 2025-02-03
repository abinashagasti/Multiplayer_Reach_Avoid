classdef Environment < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        % N = 3 % dimension of motion space
        n % number of pursuers
        m % number of evaders
        % agent positions
        pursuer_positions % n*3 
        evader_positions % m*3
        % agent_velocities 
        pursuer_speeds% n*1
        evader_speeds% m*1
        alpha % speed ratio matrix m*n
        B % barrier_matrix
        a % payoff_matrix
        V % value_matrix
        x % assignment
        evader_status % indicator to stop ode45 and restart
        % 0 - free evader
        % 1 - captured evader
        % -1 - reached target away from pursuer
    end

    methods
        
        function env = Environment(n,m,pursuer_positions,evader_positions,pursuer_speeds,evader_speeds)
            % env.N = N;
            env.n = n;
            env.m = m;
            env.pursuer_positions = pursuer_positions;
            env.evader_positions = evader_positions;
            env.pursuer_speeds = pursuer_speeds;
            env.evader_speeds = evader_speeds;
            env.alpha = env.evader_speeds*(1./env.pursuer_speeds)';
            env.B = zeros(m,n);
            env.a = zeros(m,n);
            env.V = zeros(m,n);
            env.x = zeros(m,n);
            env.evader_status = zeros(m,1);
        end

        function compute_barrier(env)
            env.B = vecnorm(env.evader_positions,2,2).^2-(vecnorm(env.pursuer_positions,2,2).^2)'.*env.alpha.^2;
        end

        function compute_value_matrix(env)
            % Takes in the position a pursuer and an evader and their speed ratios
            % Returns a mxn matrix with each entry containing V_ij

             % Compute the squared norms of each row
            norm_pursuer_squared = sum(env.pursuer_positions.^2, 2); % n x 1
            norm_evader_squared = sum(env.evader_positions.^2, 2);   % m x 1
            
            % Compute pairwise distances between pursuer_positions and evader_positions
            distances1 = sqrt(sum((reshape(env.pursuer_positions, size(env.pursuer_positions, 1), 1, size(env.pursuer_positions, 2)) - ...
                     reshape(env.evader_positions, 1, size(env.evader_positions, 1), size(env.evader_positions, 2))).^2, 3)); % n x m
            
            % Compute the matrix in case alpha=1
            value1 = 0.5 * (norm_evader_squared' - norm_pursuer_squared) ./ distances1;
            value1 = value1';
            
            % Compute pairwise differences and norms
            differences = reshape(env.evader_positions, size(env.evader_positions, 1), 1, size(env.evader_positions, 2)) - ...
                          reshape(env.pursuer_positions, 1, size(env.pursuer_positions, 1), size(env.pursuer_positions, 2));
            
            % Pairwise norm between pursuer_positions and evader_positions
            distances2 = sqrt(sum(differences.^2, 3)); % m x n
            
            % Reshape alpha for broadcasting
            alpha_squared = env.alpha.^2;
            alpha_factor = env.alpha ./ (1 - alpha_squared); % m x n
            denom_factor = 1 - alpha_squared; % m x n
            
            % Compute the first term
            first_term = sqrt(sum(((reshape(env.evader_positions, size(env.evader_positions, 1), 1, size(env.evader_positions, 2)) - ...
                                  (reshape(alpha_squared, size(alpha_squared, 1), size(alpha_squared, 2), 1) .* ...
                                  reshape(env.pursuer_positions, 1, size(env.pursuer_positions, 1), size(env.pursuer_positions, 2)))).^2), 3)) ./ denom_factor;
            
            % Compute the second term
            second_term = alpha_factor .* distances2;
            
            % Resultant matrix
            value2 = first_term - second_term;
            
            value3 = -sqrt(norm_pursuer_squared)'+sqrt(norm_evader_squared)./env.alpha;
            
            idx1 = (env.B>=0).*(env.alpha==1);
            idx2 = (env.B>=0).*(env.alpha<1);
            idx3 = (env.B<0).*(env.alpha<=1);
            v = idx1.*value1 + idx2.*value2 + idx3.*value3;
            env.V = v;
            
            % The below operation is happening in the above code in a
            % vectorized manner.
            %V=zeros(m,n);
            % for j=1:n
            %     for i=1:m
            %         if alpha(i,j)<1
            %             xc=(evader_positions(i,:)-alpha(i,j)^2*pursuer_positions(j,:))/(1-alpha(i,j)^2);
            %             rc=(alpha(i,j)/(1-alpha(i,j)^2))*norm(pursuer_positions(j,:)-evader_positions(i,:));
            %             V(i,j)=norm(xc)-rc;
            %         elseif alpha(i,j)==1
            %             V(i,j)=0.5*((norm(evader_positions(i,:))^2-norm(pursuer_positions(j,:))^2)/norm(pursuer_positions(j,:)-evader_positions(i,:)));
            %         else
            %             V(i,j)=-norm(pursuer_positions(j,:))+norm(evader_positions(i,:))/alpha(i,j);
            %         end
            %     end
            % end

            % Next we compute the payoff matrix from the value matrix.
            % env.compute_barrier();
            idx = (env.B>=0).*(env.alpha<=1);
            a = idx.*v;
            L = sum(max(a,[],2));
            a(a==0) = -L-1;
            env.a = a;
        end

        function optimal_assignment(env)
            f = env.a';
            f = f(:);
            b = ones(env.n+env.m,1);
            A = zeros(env.n+env.m,env.n*env.m);
            row_indices = (env.n+1:env.n+env.m)';
            col_start_indices = (row_indices-env.n-1)*env.n+1;
            col_indices = (0:env.n-1)+col_start_indices;
            linear_indices = sub2ind(size(A),repmat(row_indices,1,env.n),col_indices);
            A(linear_indices)=1;
            A(1:env.n,:)=repmat(eye(env.n),1,env.m);
            if env.n>=env.m
                x=linprog(-f,A(1:env.n,:),b(1:env.n),A(env.n+1:env.n+env.m,:),b(env.n+1:env.n+env.m),zeros(env.m*env.n,1));
                x=reshape(x,[env.n,env.m])';
            end
            if env.n<env.m
                x=linprog(-f,A(env.n+1:env.n+env.m,:),b(env.n+1:env.n+env.m),A(1:env.n,:),b(1:env.n),zeros(env.m*env.n,1));
                x=reshape(x,[env.n,env.m])';
            end
            env.x = x;
        end

        function [win,check_cond] = check_win(env)
            % env.compute_value_matrix();
            % env.optimal_assignment();
            check=env.a.*env.x;
            if min(check(:))<0
                win=1;
            else
                win=0;
            end
            % evader_contributions=sum(env.a.*env.x,2);
            % % Logical condition for evader_contributions <= 0
            % evader_condition = evader_contributions <= 0; % m x 1 logical vector
            % 
            % % Find indices where x(i, :) == 1 for all rows simultaneously
            % [rows, cols] = find(env.x == 1); % Get row and column indices of ones in x
            % 
            % % Create a logical mask for the specific condition
            % mask = (evader_condition(rows)).*(env.alpha(sub2ind(size(env.alpha), rows, cols))>1);
            % 
            % % Determine check_cond
            % if any(mask)
            %     check_cond = 0;
            % else
            %     check_cond = 1;
            % end

            mask = (env.alpha>1).*env.x;
            
            if any(mask(:))
                check_cond=1;
            else
                check_cond=0;
            end

            % The above vectorised code to obtain check_cond achieves the
            % below code. Here we are seeing if there is an assignment that
            % results in B_ij>=0 and alpha_ij>1. That is not something we
            % desire. So here we only check if there is a situation like
            % above. If there is we don't do the simulation. a_ij takes
            % negative values if B_ij<0 or alpha_ij>1. If it is negative
            % and B_ij<0 with alpha_ij<=1 then fine. But if the case that
            % alpha_ij>1 comes then we need to stop. 
            % check_cond = 0;
            % for i=1:length(evader_contributions)
            %     if evader_contributions(i)<=0
            %         index=find(env.x(i,:)==1);
            %         if env.alpha(i,index)<=1
            %             check_cond=0;
            %         else
            %             check_cond=1;
            %         end
            %     end
            % end
        end

        function dydt = agent_dynamics_speed(env,t,y)
        % Agent dynamics for the reach-avoid game

            i=y(7);j=y(8);
            xp=y(4:6,1);xe=y(1:3,1);
            ve=env.evader_speeds(i);
            vp=env.pursuer_speeds(j);
            if env.alpha(i,j)<1
                xc=(xe-env.alpha(i,j)^2*xp)/(1-env.alpha(i,j)^2);
                rc=(env.alpha(i,j)/(1-env.alpha(i,j)^2))*norm(xp-xe);
            end
            if env.alpha(i,j)>1
                xc=(env.alpha(i,j)^2*xp-xe)/(env.alpha(i,j)^2-1);
                rc=(env.alpha(i,j)/(env.alpha(i,j)^2-1))*norm(xp-xe);
            end
            Rc=norm(xc);
            Rp=norm(xp);
            Re=norm(xe);
            
            if env.B(i,j)>=0
                if abs(norm(xp-xe))>1e-2
                    gradV=(1/(1-env.alpha(i,j)^2))*[xc/Rc-(env.alpha(i,j)^2/(1-env.alpha(i,j)^2))*((xe-xp)/rc);
                        env.alpha(i,j)^2*(-xc/Rc+(1/(1-env.alpha(i,j)^2))*((xe-xp)/rc))];%finding the 
                    % gradient of the Value function
                    rhoe=norm(gradV(1:3,1));
                    rhop=norm(gradV(4:6,1));    
                    dydt=[-(ve/rhoe)*gradV(1:3,1);(vp/rhop)*gradV(4:6,1);0;0];
        %             dydt=[-(ve/Re)*xe;(vp/rhop)*gradV(4:6,1);0;0;0];
                else
                    dydt=zeros(8,1);
                end
            elseif env.B(i,j)<0
                if abs(norm(xe))>1e-3
                    gradV=[-(1/env.alpha(i,j))*(xe/Re);xp/Rp];
                    rhoe=norm(gradV(1:3,1));
                    rhop=norm(gradV(4:6,1)); 
        %             dydt=[(alpha/rhoe)*gradV(1:3,1);-(1/rhop)*gradV(4:6,1);0;0;0];
                    dydt=[(ve/rhoe)*gradV(1:3,1);-(vp/rhop)*gradV(4:6,1);0;0];
                else
                    dydt=zeros(8,1);
                end
            end

        end

        function update_position

        function dydt = multi_agent_dynamics(env,t,y)
            
        end

        function plot_trajectories(env)

        end

        function plot_trajectories_classic(env)
            figure
            
            plot3(0,0,0,'go',LineWidth=3,MarkerSize=8)
            hold on
            plot3(env.pursuer_positions(:,1),env.pursuer_positions(:,2),env.pursuer_positions(:,3), 'ro',LineWidth=3,MarkerSize=8)
            plot3(env.evader_positions(:,1),env.evader_positions(:,2),env.evader_positions(:,3), 'bo',LineWidth=3,MarkerSize=8)
            grid on
            xlabel('x');ylabel('y');zlabel('z');
            
            env.compute_barrier();
            env.compute_value_matrix();
            env.optimal_assignment();
            [win,check_cond] = env.check_win();

            if win==0
                disp("The pursuing team wins!")
            else
                disp("The evading team wins!")
            end
            
            % if check_cond==0
            %     T=15;
            %     z=zeros(T*100+1,9,env.m);
            %     V=zeros(M,1);
            %     for i=1:M
            %         j=find(x(i,:)==1);
            %         V(i)=Value(xp(j,:),xe(i,:),alpha(i,j));
            %         [t,y]=ode45(@agent_dynamics_speed,0:0.01:T,[xe(i,:)';xp(j,:)';env.evader_speeds(i);env.pursuer_speeds(j);env.B(i,j)]);
            %         z(:,:,i)=y;
            %         plot3(z(:,1,i),z(:,2,i),z(:,3,i),'b',LineWidth=2)
            %         plot3(z(:,4,i),z(:,5,i),z(:,6,i),'r',LineWidth=2)
            %     end 
            %     netValue=sum(V)
            % end

            T = 15;
            dt = 0.01; % Time step
            time_vector = 0:dt:T; % Time vector
            z = zeros(length(time_vector), 8, env.m); % Preallocate results matrix
            
            % Identify all active columns in x (where x(i, :) == 1)
            [row_indices, col_indices] = find(env.x == 1); % Find (i, j) pairs
            
            % Iterate over each active agent pair
            if check_cond==0
                for idx = 1:length(row_indices)
                    i = row_indices(idx);
                    j = col_indices(idx);
                
                    % Solve ODE for the current pair
                    [~, y] = ode45(@env.agent_dynamics_speed, time_vector, [env.evader_positions(i, :)';...
                        env.pursuer_positions(j, :)'; i; j]);
                    z(:, :, i) = y;
                
                    % Plot the results for this pair
                    plot3(z(:, 1, i), z(:, 2, i), z(:, 3, i), 'b', LineWidth=2);
                    hold on;
                    plot3(z(:, 4, i), z(:, 5, i), z(:, 6, i), 'r', LineWidth=2);
                end
                hold off;
            end

        end

    end
end