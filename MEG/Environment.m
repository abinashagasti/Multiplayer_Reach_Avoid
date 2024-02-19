classdef Environment
    properties
        motion_space_dimension % dimension of motion space
        evader_numbers % number of evaders in the game
        pursuer_speed % pursuer speed
        evader_speeds % array of evader speeds
        target_position % position of the target
    end
    methods
        function [env, pursuer, evaders] = Environment(N, n, v, u)
            % Constructor, assigns values and creates obstacles (optional)
            env.evader_numbers = n;
            env.motion_space_dimension = N;
            env.pursuer_speed = v;
            env.evader_speeds = u;
            env.target_position = rand(env.motion_space_dimension,1);

            % Initialize pursuer and evader object instances 
            pursuer = Pursuer(rand(env.motion_space_dimension,1), env.pursuer_speed);
            evaders = Evader.empty(env.evader_numbers,0);
            for i=1:env.evader_numbers
                evaders(i) = Evader(rand(env.motion_space_dimension,1), env.evader_speeds(i));
            end
        end
        
        function reset(env, pursuer, evaders)
            pursuer.updatePos(rand(env.motion_space_dimension,1));
            for i=1:env.evader_numbers
                evaders(i).updatePos(Evader(rand(env.motion_space_dimension,1), env.evader_speeds(i)));
            end
        end
        
        function barrier_value = barrier(env, pursuer, evaders)
            barriers = zeros(env.evader_numbers,1);
            for i=1:env.evader_numbers
                barriers(i) = norm(env.target_position - pursuer.position) - norm(env.target_position - evaders(i).position);
            end
            barrier_value = min(barriers);
        end

        function check_initialization(env, pursuer, evaders, display_info)
            if barrier(env, pursuer, evaders)>=0
                disp("win = 1")
                if display_info
                    disp("The current initialization results in evaders winning. Please reset the environment and try again.")
                end
            else
                disp("win = 0")
                if display_info
                    disp("The current initialization results in pursuers winning. You may proceed.")
                end
            end
        end

    end
end

% Helper function to check for collision with a circular obstacle
