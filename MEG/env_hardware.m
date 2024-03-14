classdef env_hardware < Environment
    properties
        pos_vrpn
        ori_vrpn
    end

    methods 

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

    end
end