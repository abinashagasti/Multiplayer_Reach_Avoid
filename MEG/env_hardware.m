classdef env_hardware < Environment
    properties
        pos_vrpn
        ori_vrpn
    end

    methods 
        function env = env_hardware(N,n,v,u,r,pursuer_poisition,evader_positions, target_position)
            env= env@Environment(N,n,v,u,r,pursuer_poisition,evader_positions,target_position);
            rosinit
            robotpos_p = rossubscriber("/vrpn_client_node/pursuer/pose",@env.pursuer.callback,"DataFormat","struct");
            robotpos_e = Subscriber.empty(env.evader_numbers,0);
            for i = 1:env.envader_numbers
                evader = env.evaders(i);
                robotpos_e(i) = rossubscriber("/vrpn_client_node/"+env.evaders(i).name+"/pose",@evader.callback,"DataFormat","struct");
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

    end
end