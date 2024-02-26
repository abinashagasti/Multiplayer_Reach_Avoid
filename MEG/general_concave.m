clear
clc

syms xt yt xp yp xi yi psi theta r

m = (r*(cos(theta)-cos(psi))+xp-xi)/(yi-yp+r*(sin(psi)-sin(theta)));

c = (yi^2-xi^2+xp^2-yp^2+2*r*yi*sin(psi)-2*r*yp*sin(theta)+2*r*xi*cos(psi)-2*r*xp*cos(theta))/(2*(yi-yp+r*(sin(psi)-sin(theta))));

%%

denom = simplify(1+m^2,1000);
denom = sqrt(denom);

d = simplify(m*xt-yt+c,1000)
d = simplify(d/denom, 1000);

%% 

p_pos = [0,1]';
e_pos = [3,3]';

xc = 1.5;
yc = 2;
m = 2/3;

x_intercept = (m*(yc)+xc)/(1+m^2);
y_intercept = m*(x_intercept);

dir = [x_intercept - 3, y_intercept - 3];
dir = (1/norm(dir,2))*(dir);
psi = atan2(dir(2),dir(1));

%%

xe = 3; ye = 3;
xp = 0; yp = 1;
r = 0.01;

theta_values = linspace(-pi+0.001,pi,100);
for i=1:100
    theta = theta_values(i);
    c(i) = (ye^2-xe^2+xp^2-yp^2+2*r*ye*sin(psi)-2*r*yp*sin(theta)+2*r*xe*cos(psi)-2*r*xp*cos(theta))/(2*(ye-yp+r*(sin(psi)-sin(theta))));
    d(i) = c(i)/(1+m^2);
end
plot(theta_values,d)

%% 

dir = [x_intercept - xp, y_intercept - yp];
dir = (1/norm(dir,2))*(dir);
theta_star = atan2(dir(2),dir(1))
