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
evader_positions = [5,3,2;2,4,4];
% evader_positions = evader_positions(:,1);
target_position = [0,0];

[env, pursuer, evaders] = Environment(N,n,v,u,r,pursuer_position,evader_positions,target_position); % Initialize the environment

win = env.check_initialization(pursuer, evaders, false); % Check win scenario


%%

env.plot_trajectories(pursuer, evaders)
