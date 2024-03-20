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

% pursuer_position = [-1.8311;-1.2004];
% evader_positions = [-2.4013,-0.6859,8;3.007,4.1065,-5];

% pursuer_position = [0.5;0.5]
% evader_positions = [4,4;3,5]

pursuer_position = -4+8*rand(N,1)
evader_positions = -10+20*rand(N,n)
evader_positions = evader_positions(:,1:n);
target_position = [0;0];

env = Environment(N,n,v,u,r,pursuer_position,evader_positions,target_position); % Initialize the environment
% env = env_hardware(N,n,v,u,r,pursuer_position,evader_positions,target_position);
win = env.check_initialization(env.evaders, false); % Check win scenario
env.plot_current_positions()

%%

if win
    % set(1,'DefaultFigureWindowStyle','docked')
    [win, pursuer_position_traj, evader_positions_traj] = env.obtain_trajectories('standard');
    % env.plot_trajectories(pursuer_position_traj, evader_positions_traj)
    % xlim([-1,8]);
    % ylim([-1e-16,1e-16])
else 
    % disp('reinitialize')
    main
end

%%
