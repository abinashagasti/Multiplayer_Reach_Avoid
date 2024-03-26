classdef Environment < handle
    properties
        motion_space_dimension % dimension of motion space
        evader_numbers % number of evaders in the game
        pursuer_speed % pursuer speed
        evader_speeds % array of evader speeds
        target_position % position of the target
        timestep % timestep value
        captured_evaders % boolean array containing capture status of evaders
        pursuer % pursuer object
        evaders % evader object array
        capture_tolerance % tolerance for point capture
        alpha % speed ratio array
    end
    methods

        function env = Environment(N, n, v, u, r, pursuer_position, evader_positions, target_position)
            env.evader_numbers = n;
            env.motion_space_dimension = N;
            env.pursuer_speed = v;
            env.evader_speeds = u;
            env.timestep = r;
            env.captured_evaders = boolean(zeros(1,env.evader_numbers));
            env.evaders = Evader.empty(env.evader_numbers,0);
            env.capture_tolerance = 0.05;
            env.alpha = zeros(1,n);
            for i=1:n
                env.alpha(i)=env.evader_speeds(i)/env.pursuer_speed;
            end

            % Initialize pursuer and evader object instances 
            if ~isempty(pursuer_position)
                env.pursuer = Pursuer(pursuer_position, env.pursuer_speed);
            else
                env.pursuer = Pursuer(rand(env.motion_space_dimension,1), env.pursuer_speed);
            end
            
            if ~isempty(evader_positions)
                for i=1:env.evader_numbers
                   env.evaders(i) = Evader(evader_positions(:,i), env.evader_speeds(i), i);
                end
            else
                for i=1:env.evader_numbers
                   env.evaders(i) = Evader(rand(env.motion_space_dimension,1), env.evader_speeds(i), i);
                end
            end
            if ~isempty(target_position)
                env.target_position = target_position;
            else
                env.target_position = rand(env.motion_space_dimension,1);
            end
        end

        function update_target(env, target_position)
            sz_input = size(target_position);
            sz_output = size(env.target_position);
            if sz_input~=sz_output
                if sz_input==fliplr(sz_output)
                    env.target_position = target_position';
                else
                    error('Wrong position dimensions.')
                end
            else
                env.target_position = target_position;
            end
        end
        
        function reset(env)
            env.pursuer.updatePos(rand(env.motion_space_dimension,1));
            for i=1:env.evader_numbers
                env.evaders(i).updatePos(Evader(rand(env.motion_space_dimension,1), env.evader_speeds(i)));
            end
        end

        function update_Environment(env,pursuer_position, evader_positions, target_position)
            env.pursuer.updatePos(pursuer_position)
            for i=1:env.evader_numbers
                env.evaders(i).updatePos(evader_positions(:,i));
            end
            env.target_position = target_position;
        end
        
        function barrier_value = barrier(env, evaders)
            barriers = zeros(length(evaders),1);
            % disp(length(evaders))
            for i=1:length(evaders)
                barriers(i) = norm(env.target_position - evaders(i).position)^2 - env.alpha(i)^2*norm(env.target_position - env.pursuer.position)^2;
            end
            barrier_value = min(barriers);
        end

        function win = check_initialization(env,evaders,display_info)
            if env.barrier(evaders)<0
                win = false;
                if display_info
                    disp("The current initialization results in evaders winning. Please reset the environment and try again.")
                end
            else
                win = true;
                if display_info
                    disp("The current initialization results in pursuers winning. You may proceed.")
                end
            end
        end

        function evader_names = return_evader_names(env)
            evader_names = strings(1,env.evader_numbers);
            for i=1:env.evader_numbers
                env.evaders(i).name = "evader"+int2str(env.evaders(i).index);
                evader_names(i) = env.evaders(i).name;
            end
        end

        function plot_current_positions(env)
            hold on
            for i=1:env.evader_numbers
                plot(env.evaders(i).position(1), env.evaders(i).position(2), '.', 'color', 'b', 'MarkerSize', 30)
            end
            plot(env.target_position(1), env.target_position(2), '.', 'color', 'g', 'MarkerSize', 30)
            plot(env.pursuer.position(1), env.pursuer.position(2), '.', 'color', 'r', 'MarkerSize', 30)
            hold off
        end

        function evader_positions = return_evader_positions(env,evaders)
            evader_positions = zeros(env.motion_space_dimension,length(evaders));
            for i=1:length(evaders)
                evader_positions(:,i) = evaders(i).position;
            end
        end

        function evader_velocities = return_evader_velocities(env)
            evader_velocities = zeros(env.motion_space_dimension,env.evader_numbers);
            for i=1:env.evader_numbers
                win = env.check_initialization(env.evaders(i),false);
                % win = true;
                evader_velocities(:,i) = env.evaders(i).heading_velocity(env.pursuer.position, env.target_position,win,env.alpha(i));
            end
        end

        function updateTermination(env)
            for i=1:env.evader_numbers
                if norm(env.pursuer.position-env.evaders(i).position)<env.capture_tolerance
                    env.captured_evaders(i) = true;
                end
            end
        end

        function done = step(env, objective_function)
            env.updateTermination();
            % if any(env.captured_evaders)
            %     done = true;
            %     return
            % end
            done = false;
            if all(env.captured_evaders)
                done = true;
                return
            end
            evader_list = 1:env.evader_numbers;
            if ~done
                for i=evader_list(~env.captured_evaders)
                    if norm(env.evaders(i).position- env.target_position)<env.capture_tolerance && norm(env.pursuer.position-env.target_position)>2*env.capture_tolerance
                        done = true;
                        disp_message = strcat(env.evaders(i).name,' reached the target.');
                        disp(disp_message)
                        return
                    end
                end
            end
            win = env.check_initialization(env.evaders(evader_list(~env.captured_evaders)),false);
            pursuer_velocity = env.pursuer.heading_velocity(env.return_evader_positions(env.evaders(evader_list(~env.captured_evaders))), env.target_position, env.timestep, win, objective_function);
            evader_velocities = env.return_evader_velocities();
            env.pursuer.updatePos(env.pursuer.position + env.timestep*pursuer_velocity);
            for i=evader_list(~env.captured_evaders)
                env.evaders(i).updatePos(env.evaders(i).position + env.timestep*evader_velocities(:,i));
            end
        end


        function [win, pursuer_positions_traj, evader_positions_traj] = obtain_trajectories(env, objective_function)
            env.plot_current_positions()
            hold on
            done = false; t = 1;
            while ~done
                pursuer_positions_traj(:,t) = env.pursuer.position;
                for i=1:env.evader_numbers
                    evader_positions_traj(:,t,i) = env.evaders(i).position;
                end
                if mod(t,10)
                    pause(0.01)
                    plot(pursuer_positions_traj(1,:),pursuer_positions_traj(2,:),'r')
                    hold on
                    for i=1:env.evader_numbers
                        plot(evader_positions_traj(1,:,i),evader_positions_traj(2,:,i),'b')
                    end
                    env.plot_current_positions()
                end
                done = env.step(objective_function);
                t = t+1;
                % win = env.check_initialization(env.evaders,false);
                % if win==0
                %     done = true;
                %     return
                % end
            end
            pursuer_positions_traj(:,t) = env.pursuer.position;
            for i=1:env.evader_numbers
                evader_positions_traj(:,t,i) = env.evaders(i).position;
            end
            if all(env.captured_evaders)
                win=true;
            else
                win=false;
            end
            plot(pursuer_positions_traj(1,:),pursuer_positions_traj(2,:),'r')
            hold on
            for i=1:env.evader_numbers
                plot(evader_positions_traj(1,:,i),evader_positions_traj(2,:,i),'b')
            end
            env.plot_current_positions()
        end

        function plot_trajectories(env, pursuer_positions_traj, evader_positions_traj)
            % hold on
            plot(pursuer_positions_traj(1,:),pursuer_positions_traj(2,:),'r')
            for i=1:env.evader_numbers
                plot(evader_positions_traj(1,:,i),evader_positions_traj(2,:,i),'b')
            end
            env.plot_current_positions()
        end


    end
end

% Helper function to check for collision with a circular obstacle
