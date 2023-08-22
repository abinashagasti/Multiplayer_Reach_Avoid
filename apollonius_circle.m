%% Code to find apollonius circle parameters given xp and xe


clear
clc

xe=[-1,1];
xp=[1,1];
alpha=0.5;


if alpha<1
    xc=(xe-alpha^2.*xp)/(1-alpha^2)
    rc=(alpha/(1-alpha^2))*sqrt((xp(1)-xe(1))^2+(xp(2)-xe(2))^2)
    
    B=(xe(1)^2+xe(2)^2)-alpha^2*(xp(1)^2+xp(2)^2)
    
    theta=7*pi/4;
    xa=xc+rc.*[cos(theta),sin(theta)]
else
    xc=(xe-alpha^2.*xp)/(1-alpha^2)
    rc=-(alpha/(1-alpha^2))*sqrt((xp(1)-xe(1))^2+(xp(2)-xe(2))^2)
    
    B=(xe(1)^2+xe(2)^2)-alpha^2*(xp(1)^2+xp(2)^2)
    
    theta=7*pi/4;
    xa=xc+rc.*[cos(theta),sin(theta)]
end