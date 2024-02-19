%% 
clear
clc
close all

%% Initialize the environment and the agents

N = 2; % Dimension of motion space (2D plane)
n = 1; % Number of evaders
v = 1; % Pursuer speed
u = ones(1,n); % Evader speeds

[env, pursuer, evaders] = Environment(N,n,v,u); % Initialize the environment

%%

env.check_initialization(pursuer, evaders, false)