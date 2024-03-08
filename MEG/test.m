%% Testing usage of convex minimization

clear
clc

y = 5;

xmin = -10; xmax = 10;

A = [1;-1]; b = [xmax;-xmin];

fun = @(x) test_fun(x,y);

fmincon(fun,0,A,b)

fun(3)

%% Testing pursuer class functions

clear
clc

pursuer = Pursuer([1;0],1);
evader_positions = [5,3,-1;2,4,6];
target_position = [0;0];

psis = zeros(3,1);
thetas = zeros(3,1);
for i=1:3
    [psis(i), thetas(i)] = pursuer.optimal_headings_Ei(evader_positions(:,i),target_position);
end

% psis
% thetas

[theta_min, theta_max] = pursuer.concave_domain(evader_positions, target_position);

%% 

e