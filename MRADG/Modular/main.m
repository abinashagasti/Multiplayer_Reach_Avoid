%%

clear
clc
close all

%%

% dimension of motion space = 3 
n=randi([3 4],1); % number of pursuers
m=randi([n-2 n],1); % number of evaders

% randomized initializations of pursuer and evader positions
pursuer_positions = -10+20*rand(n,3); 
evader_positions = -15+30*rand(m,3);

% randomized initializations of pursuer and evader velocities
pursuer_speeds = 1.5+rand(n,1);
evader_speeds = 1+rand(m,1);

env = Environment(n,m,pursuer_positions,evader_positions,pursuer_speeds,evader_speeds);

%% 
env.plot_trajectories()
