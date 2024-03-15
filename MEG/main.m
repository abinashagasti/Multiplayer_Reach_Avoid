%% 
clear
clc
close all

%% Initialize the environment and the agents

N = 2; % Dimension of motion space (2D plane)
n = 1; % Number of evaders
v = 1; % Pursuer speed
u = ones(1,n); % Evader speeds
r = 0.01; % Timestep

pursuer_position = [1;0];
evader_positions = [3,4,2;-0,0,4];
evader_positions = evader_positions(:,1);
target_position = [0,0];

% env = Environment(N,n,v,u,r,pursuer_position,evader_positions,target_position); % Initialize the environment
env = env_hardware(N,n,v,u,r,pursuer_position,evader_positions,target_position);
% win = env.check_initialization(false); % Check win scenario


%%

% env.plot_trajectories()
% xlim([-1,8]);
% ylim([-1e-16,1e-16])