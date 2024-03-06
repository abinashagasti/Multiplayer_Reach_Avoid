%% 
clear
clc
close all

%% Initialize the environment and the agents

N = 2; % Dimension of motion space (2D plane)
n = 3; % Number of evaders
v = 1; % Pursuer speed
u = ones(1,n); % Evader speeds

[env, pursuer, evaders] = Environment(N,n,v,u); % Initialize the environment

win = env.check_initialization(pursuer, evaders, false); % Check win scenario

%% Setting some positions 

pursuer.updatePos([1,0]);
evaders(1).updatePos([5,2]);
evaders(2).updatePos([3,4]);
evaders(3).updatePos([-1,6]);
env.update_target([4,3]);

%%

evader_positions = env.return_evader_positions(evaders);

[theta_min, theta_max] = pursuer.concave_domain(evader_positions, env.target_position)