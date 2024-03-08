classdef Environment < handle
    properties
        motion_space_dimension % dimension of motion space
        evader_numbers % number of evaders in the game
        pursuer_speed % pursuer speed
        evader_speeds % array of evader speeds
        target_position % position of the target
        timestep % timestep value
        captured_evaders % boolean array containing capture status of evaders
    end
    methods

        function [env, pursuer, evaders] = Environment(N, n, v, u, r, pursuer_position, evader_positions, target_position)
            env.evader_numbers = n;
            env.motion_space_dimension = N;
            env.pursuer_speed = v;
            env.evader_speeds = u;
            env.timestep = r;
            env.captured_evaders = boolean(zeros(1,env.evader_numbers));

            % Initialize pursuer and evader object instances 
            if ~isempty(pursuer_position)
                pursuer = Pursuer(pursuer_position, env.pursuer_speed);
            else
                pursuer = Pursuer(rand(env.motion_space_dimension,1), env.pursuer_speed);
            end
            evaders = Evader.empty(env.evader_numbers,0);
            if ~isempty(evader_positions)
                for i=1:env.evader_numbers
                   evaders(i) = Evader(evader_positions(:,i), env.evader_speeds(i), i);
                end
            else
                for i=1:env.evader_numbers
                   evaders(i) = Evader(rand(env.motion_space_dimension,1), env.evader_speeds(i), i);
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
        
        function reset(env, pursuer, evaders)
            pursuer.updatePos(rand(env.motion_space_dimension,1));
            for i=1:env.evader_numbers
                evaders(i).updatePos(Evader(rand(env.motion_space_dimension,1), env.evader_speeds(i)));
            end
        end

        function update_Environment(env, pursuer, evaders, pursuer_position, evader_positions, target_position)
            pursuer.updatePos(pursuer_position)
            for i=1:env.evader_numbers
                evaders(i).updatePos(evader_positions(:,i));
            end
            env.target_position = target_position;
        end
        
        function barrier_value = barrier(env, pursuer, evaders)
            barriers = zeros(env.evader_numbers,1);
            for i=1:env.evader_numbers
                barriers(i) = norm(env.target_position - evaders(i).position) - norm(env.target_position - pursuer.position);
            end
            barrier_value = min(barriers);
        end

        function win = check_initialization(env, pursuer, evaders, display_info)
            if barrier(env, pursuer, evaders)<0
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

        function plot_current_positions(env, pursuer, evaders)
            plot(pursuer.position(1), pursuer.position(2), '.', 'color', 'r', 'MarkerSize', 30)
            hold on
            for i=1:env.evader_numbers
                plot(evaders(i).position(1), evaders(i).position(2), '.', 'color', 'b', 'MarkerSize', 30)
            end
            plot(env.target_position(1), env.target_position(2), '.', 'color', 'g', 'MarkerSize', 30)
            hold off
        end

        function evader_positions = return_evader_positions(env, evaders)
            evader_positions = zeros(env.motion_space_dimension,env.evader_numbers);
            for i=1:env.evader_numbers
                evader_positions(:,i) = evaders(i).position;
            end
        end

        function evader_velocities = return_evader_velocities(env, pursuer, evaders)
            evader_velocities = zeros(env.motion_space_dimension,env.evader_numbers);
            for i=1:env.evader_numbers
                evader_velocities(:,i) = evaders(i).heading_velocity(pursuer.position, env.target_position);
            end
        end

        function updateTermination(env, pursuer, evaders)
            env.captured_evaders = boolean(zeros(1,env.evader_numbers));
            for i=1:env.evader_numbers
                if norm(pursuer.position-evaders(i).position)<0.1
                    env.captured_evaders(i) = true;
                end
            end
        end

        function done = step(env, pursuer, evaders)
            env.updateTermination(pursuer,evaders);
            if ~any(~env.captured_evaders)
                done = true;
                return
            end
            evader_list = 1:env.evader_numbers;
            pursuer_velocity = pursuer.heading_velocity(env.return_evader_positions(evaders(evader_list(~env.captured_evaders))), env.target_position, env.timestep);
            evader_positions = env.return_evader_positions(evaders);
            evader_velocities = env.return_evader_velocities(pursuer, evaders);
            pursuer.updatePos(pursuer.position + env.timestep*pursuer_velocity);
            for i=1:evader_list(~env.captured_evaders)
                evaders(i).updatePos(evader_positions(:,i) + env.timestep*evader_velocities(:,i));
            end
            done = false;
        end


        function plot_trajectories(env, pursuer, evaders)
            env.plot_current_positions(pursuer, evaders)
            hold on
            done = false; t = 1;
            while ~done
                pursuer_positions_traj(:,t) = pursuer.position;
                for i=1:env.evader_numbers
                    evader_positions_traj(:,t,i) = evaders(i).position;
                end
                done = env.step(pursuer,evaders);
                t = t+1;
            end
            pursuer_positions_traj(:,t) = pursuer.position;
            for i=1:env.evader_numbers
                evader_positions_traj(:,t,i) = evaders(i).position;
            end
            plot(pursuer_positions_traj(1,:),pursuer_positions_traj(2,:),'r')
            for i=1:env.evader_numbers
                plot(evader_positions_traj(1,:,i),plot(evader_positions_traj(2,:,i)),'b')
            end
            env.plot_current_positions(pursuer, evaders)
        end


    end
end

% Helper function to check for collision with a circular obstacle
