%% 
clear
clc
close all

%% Initialize the environment and the agents

N = 2; % Dimension of motion space (2D plane)
n = 3; % Number of evaders
v = 1; % Pursuer speed
u = ones(1,n); % Evader speeds
r = 0.01; % Timestep

pursuer_position = [1;0];
evader_positions = [5,3,-1;2,4,6];
target_position = [0,0];

[env, pursuer, evaders] = Environment(N,n,v,u,r,pursuer_position,evader_positions,target_position); % Initialize the environment

win = env.check_initialization(pursuer, evaders, false); % Check win scenario


%%

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

%%

[thetamin, thetamax] = pursuer.concave_domain(evader_positions, target_position)