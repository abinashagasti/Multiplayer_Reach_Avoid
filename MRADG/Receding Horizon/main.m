%%

clear
clc
close all
set(0, 'DefaultFigureWindowStyle', 'docked');

%

% dimension of motion space = 3 
n=randi([3 4],1); % number of pursuers
m=randi([max(1,n-2) n],1); % number of evaders

% randomized initializations of pursuer and evader positions
pursuer_positions = -10+20*rand(n,3); 
evader_positions = -15+30*rand(m,3);

% randomized initializations of pursuer and evader velocities
pursuer_speeds = 1.5+rand(n,1);
evader_speeds = 1+rand(m,1);

%%
env = Environment(n,m,pursuer_positions,evader_positions,pursuer_speeds,evader_speeds);

env.check_win()
env.win
env.plot_current_positions()
done = false; t = 1;
time_chunk = 1000;
pursuer_traj = zeros(3,time_chunk,env.n);
evader_traj= zeros(3,time_chunk,env.m);

%%
while ~done
    pursuer_traj(:,t,:) = reshape([env.pursuers.position], 3, 1, []);
    evader_traj(:,t,:) = reshape([env.evaders.position], 3, 1, []);
    % if mod(t,10)
        % pause(0.01)
        % hold on
        % for j=1:env.n
        %     plot3(pursuer_traj(1,:,j),pursuer_traj(2,:,j),pursuer_traj(3,:,j),'r')
        % end
        % for i=1:env.m
        %     plot3(evader_traj(1,:,j),evader_traj(2,:,j),evader_traj(3,:,j),'b')
        % end
        % env.plot_current_positions()
    % end
    done = env.step();
    t = t+1;
    if t>time_chunk
        pursuer_traj(:, end+time_chunk, :) = 0;
        evader_traj(:, end+time_chunk, :) = 0;
    end
end

%%

% xp = env.pursuers(1).position;
% xe = env.evaders(1).position;alpha=env.alpha(1,1);
% xc = (xe-alpha^2*xp)/(1-alpha^2);
% Rc = norm(xc);
% rc = (alpha/(1-alpha^2))*norm(xp-xe);
% norm((1-rc/Rc)*xc-xp)+norm((1-rc/Rc)*xc-xe)


%%
pursuer_traj = pursuer_traj(:, 1:t, :);
evader_traj= evader_traj(:, 1:t, :);

pursuer_traj(:,t,:) = reshape([env.pursuers.position], 3, 1, []);
evader_traj(:,t,:) = reshape([env.evaders.position], 3, 1, []);

hold on
for j=1:env.n
    plot3(pursuer_traj(1,:,j),pursuer_traj(2,:,j),pursuer_traj(3,:,j),'r')
end
for i=1:env.m
    plot3(evader_traj(1,:,i),evader_traj(2,:,i),evader_traj(3,:,i),'b')
end
% env.plot_current_positions()
