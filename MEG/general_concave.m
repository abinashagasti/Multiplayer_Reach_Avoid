%% Initialize positions
clear
clc

xp = 2; yp = 2;
xe = 1; ye = -2;
xt = 4; yt = 3;

% Find angle made by positive EP line

delta = atan2(yp-ye,xp-xe);
% delta = delta*(180/pi)

B = norm([xe,ye]-[xt,yt])-norm([xp,yp]-[xt,yt]);

xm = (xe+xp)/2; ym = (ye+yp)/2;
m = (ye-yp)/(xe-xp);

xI = (m*(ym-yt)+xm+m^2*xt)/(1+m^2);
yI = (m*(xm-xt)+yt+m^2*ym)/(1+m^2);

psi = atan2(yI-ye,xI-xe);
psi_local = psi - delta;

theta_local = pi-psi_local;
theta = atan2(yI-yp, xI-xp);
mod(theta - delta - theta_local, 2*pi) % should be equal to zero

% Here we verify that actual theta^star obtained via direction
% towards interception point and the theta_local^star obtained via local
% coordinate calculations are consistent with each other. So we can use
% the global version of the d function without any issue. 

Tx = abs(m*yt+xt-m*ym-xm)/sqrt(1+m^2);
Ty = abs(yt-m*xt-ye-m*xe)/sqrt(1+m^2);
k = norm([xp,yp]-[xe,ye])/2;
r = 0.01;
theta_values = linspace(theta-pi/2, theta+pi/2, 100);
d = zeros(length(theta_values),1);

for i=1:length(theta_values)
    angle = theta_values(i);
    d(i) = Tx + (r/(2*k))*sqrt(Ty^2+k^2)*(-1-cos(angle+psi-2*delta));
end

plot(theta_values, d)