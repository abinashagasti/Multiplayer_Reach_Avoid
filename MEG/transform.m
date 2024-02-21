function wheel_omega= transform(robot_vel,r,d)
%TRANFORM Summary of this function goes here
%   Detailed explanation goes here

jacobian = (1/r)*[-d 1 0;-d -0.5 -sin(pi/3);-d -0.5 +sin(pi/3)];
wheel_omega = jacobian*robot_vel;
    
end

