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
        pursuers % pursuer object array
        evaders % evader object array
        timestep % timestep for simulation
        tolerance % capture and invasion tolerance
        win % 0 if pursuers win and 1 otherwise
        check_cond % 1 if a pursuer-evader matching is 
        % such that alpha>1, 0 otherwise
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
            env.timestep = 0.01;
            env.tolerance = 0.01;
            env.alpha = env.evader_speeds*(1./env.pursuer_speeds)';
            env.B = zeros(m,n);
            env.a = zeros(m,n);
            env.V = zeros(m,n);
            env.x = zeros(m,n);
            env.pursuers = Pursuer.empty(env.n,0);
            env.evaders = Evader.empty(env.m,0);
            for j=1:env.n
                env.pursuers(j) = Pursuer(env.pursuer_positions(j,:)',env.pursuer_speeds(j),j);
            end
            for i=1:env.m
                env.evaders(i) = Evader(env.evader_positions(i,:)',env.evader_speeds(i),i);
            end
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
            for i=1:env.m
                if any(find(x(i,:)==1))
                    env.evaders(i).pursuer = find(x(i,:)==1);
                end
            end
            for j=1:env.n
                if any(find(x(:,j)==1))
                    env.pursuers(j).evader = find(x(:,j)==1);
                    env.pursuers(j).status = -1;
                end
            end
        end

        function check_win(env)
            env.compute_barrier()
            env.compute_value_matrix();
            env.optimal_assignment();
            [rows, columns] = find(env.x==1);
            for i=1:length(rows)
                env.evaders(rows(i)).pursuer = columns(i);
            end
            for j=1:length(columns)
                env.pursuers(columns(j)).evader = rows(j);
            end
            check=env.a.*env.x;
            if min(check(:))<0
                env.win=1;
            else
                env.win=0;
            end 

            mask = (env.alpha>1).*env.x;
            
            if any(mask(:))
                env.check_cond=1;
            else
                env.check_cond=0;
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

        function updateStatus(env)
            for i=1:env.m
                j = env.evaders(i).pursuer;
                if norm(env.pursuers(j).position-env.evaders(i).position)<env.tolerance ...
                        && norm(env.evaders(i).position)>env.tolerance
                    env.evaders(i).status = 1;
                    env.pursuers(j).status = 1;
                elseif norm(env.evaders(i).position)<env.tolerance
                    env.evaders(i).status = -1;
                end
            end
        end

        function done = terminationStatus(env)
            if ~env.win
                done = all([env.evaders.status]==1);
            else
                winning_evaders = find(sum(env.a.*env.x,2)<=0);
                if isempty(find(-[env.evaders(winning_evaders).status],1))
                    done = 0;
                else
                    done = all(find([env.evaders(winning_evaders).status]==-1));
                end
            end
        end

        function done = step(env)
            env.updateStatus();
            % if any(env.captured_evaders)
            %     done = true;
            %     return
            % end
            done = env.terminationStatus();
            for i=find([env.evaders.status]==0)
                j=find(env.x(i,:)==1);
                env.evaders(i).updatePos(env.evaders(i).position+env.timestep*...
                    env.evaders(i).return_velocity(env.pursuers(j),env.B(i,j),env.tolerance));
            end
            for j=find([env.pursuers.status]==-1)
                i=find(env.x(:,j)==1);
                env.pursuers(j).updatePos(env.pursuers(j).position+env.timestep*...
                    env.pursuers(j).return_velocity(env.evaders(i),env.B(i,j),env.tolerance));
            end
        end

        function plot_current_positions(env)
            plot3(0, 0, 0, '.', 'color', 'g', 'MarkerSize', 30)
            hold on
            for i=1:env.m
                plot3(env.evaders(i).position(1), env.evaders(i).position(2), env.evaders(i).position(3), '.', 'color', 'b', 'MarkerSize', 30)
            end
            for j=1:env.n
                plot3(env.pursuers(j).position(1), env.pursuers(j).position(2), env.pursuers(j).position(3), '.', 'color', 'r', 'MarkerSize', 30)
            end
            % xmin = min(min([env.pursuer_positions(:,1);env.evader_positions(:,1)]),0);
            % xmax = max(max([env.pursuer_positions(:,1);env.evader_positions(:,1)]),0);
            % ymin = min(min([env.pursuer_positions(:,2);env.evader_positions(:,2)]),0);
            % ymax = max(max([env.pursuer_positions(:,2);env.evader_positions(:,2)]),0);
            % zmin = min(min([env.pursuer_positions(:,3);env.evader_positions(:,3)]),0);
            % zmax = max(max([env.pursuer_positions(:,3);env.evader_positions(:,3)]),0);
            % xlim([xmin-0.2*abs(xmin), xmax+0.2*abs(xmax)]);
            % ylim([ymin-0.2*abs(ymin), ymax+0.2*abs(ymax)]);
            % zlim([zmin-0.2*abs(zmin), zmax+0.2*abs(zmax)]);
            xlabel('X-axis');
            ylabel('Y-axis');
            zlabel('Z-axis');
            grid on;
            hold off
        end

        function [pursuer_traj,evader_traj,t] = plot_trajectories(env)
            env.check_win()
            env.plot_current_positions()
            done = false; t = 1;
            time_chunk = 1000;
            pursuer_traj = zeros(3,time_chunk,env.n);
            evader_traj= zeros(3,time_chunk,env.m);
            
            while ~done
                pursuer_traj(:,t,:) = reshape([env.pursuers.position], 3, 1, env.n);
                evader_traj(:,t,:) = reshape([env.evaders.position], 3, 1, env.m);
                if mod(t,10)
                    pause(0.01)
                    hold on
                    for j=1:env.n
                        plot3(pursuer_traj(1,:,j),pursuer_traj(2,:,j),pursuer_traj(3,:,j),'r')
                    end
                    for i=1:env.m
                        plot3(evader_traj(1,:,j),evader_traj(2,:,j),evader_traj(3,:,j),'b')
                    end
                    env.plot_current_positions()
                end
                done = env.step();
                t = t+1;
                if t>time_chunk
                    pursuer_traj(:, end+time_chunk, :) = 0;
                    evader_traj(:, end+time_chunk, :) = 0;
                end
                pursuer_traj = pursuer_traj(:, 1:t, :);
                evader_traj= evader_traj(:, 1:t, :);
            end

            pursuer_traj(:,t, :) = reshape([env.pursuers.position], 3, 1, env.n);
            evader_traj(:,t,i) = reshape([env.evaders.position], 3, 1, env.m);
            
            hold on
            for j=1:env.n
                plot3(pursuer_traj(1,:,j),pursuer_traj(2,:,j),pursuer_traj(3,:,j),'r')
            end
            for i=1:env.m
                plot3(evader_traj(1,:,j),evader_traj(2,:,j),evader_traj(3,:,j),'b')
            end
            env.plot_current_positions()

        end
        

    end
end