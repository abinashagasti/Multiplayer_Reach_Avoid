%% 
clear
clc
close all

%% Initialize the environment and the agents

N = 2; % Dimension of motion space (2D plane)
n = 2; % Number of evaders
v = 1; % Pursuer speed
u = ones(1,n); % Evader speeds
r = 0.01; % Timestep

pursuer_position = [1;0];
evader_positions = [5,3;2,4];
% evader_positions = evader_positions(:,1);
target_position = [0,0];

[env, pursuer, evaders] = Environment(N,n,v,u,r,pursuer_position,evader_positions,target_position); % Initialize the environment

win = env.check_initialization(pursuer, evaders, false); % Check win scenario


%%

% env.plot_trajectories(pursuer, evaders)

env.plot_current_positions(pursuer, evaders)
hold on
done = false; t = 1;
for i=1:550
    pursuer_positions_traj(:,t) = pursuer.position;
    for i=1:n
        evader_positions_traj(:,t,i) = evaders(i).position;
    end
    done = env.step(pursuer,evaders);
    t = t+1;
end
pursuer_positions_traj(:,t) = pursuer.position;
for i=1:n
    evader_positions_traj(:,t,i) = evaders(i).position;
end
plot(pursuer_positions_traj(1,:),pursuer_positions_traj(2,:),'r')
for i=1:env.evader_numbers
    plot(evader_positions_traj(1,:,1),evader_positions_traj(2,:,1),'b')
end

%%

% env.updateTermination(pursuer,evaders);
% if ~any(~env.captured_evaders)
%     done = true;
%     return
% end
% evader_list = 1:env.evader_numbers;
% pursuer_velocity = pursuer.heading_velocity(env.return_evader_positions(evaders(evader_list(~env.captured_evaders))), env.target_position, env.timestep);
% evader_positions = env.return_evader_positions(evaders);
% evader_velocities = env.return_evader_velocities(pursuer, evaders);
% pursuer.updatePos(pursuer.position + env.timestep*pursuer_velocity);
% for i=evader_list(~env.captured_evaders)
%     evaders(i).updatePos(evaders(i).position + env.timestep*evader_velocities(:,i));
% end
% done = false;
% pursuer.position
% for i=1:env.evader_numbers
%     evaders(i).position
% end