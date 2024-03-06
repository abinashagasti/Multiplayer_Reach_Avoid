% Verification of concavity in local coordinates
clear
clc
close all 

k = 2;
xt = 6; yt = -7;
r = 0.1;

psi = atan2(yt,k);

% Optimal theta should be 2.214 rads (126.89 degrees)

theta_values = linspace(pi-psi-pi/2, pi-psi+pi/2, 200);
d = zeros(length(theta_values),1);
d_approx = d;

for i=1:length(theta_values)
    theta = theta_values(i);
    d_num = (2*k+r*(cos(theta)-cos(psi)))*xt+r*(sin(theta)-sin(psi))*yt-r*k*(cos(theta)+cos(psi));
    d_num = abs(d_num);
    d_den = 4*k^2+4*k*r*(cos(theta)-cos(psi))+r^2*(2-2*cos(theta)*cos(psi)-2*sin(theta)*sin(psi));
    d_den = sqrt(d_den);
    d(i) = d_num/d_den;
    d_approx(i) = xt+(r/(2*k))*sqrt(yt^2+k^2)*(-1-cos(theta+psi));
end

plot(theta_values, d)
hold on
plot(theta_values, d_approx)

% This code section verifies that if we consider a local coordinate system,
% i.e., P and E on the x-axis and their heading angles calculated
% anticlockwise from the x-axis, then the actual d function seems fully
% concave (not only the approximation). Further the maxima point seems to
% be what is expected, i.e., theta^star = pi - psi^star. And comparing the
% actual d function with the approximated one shows that they are almost
% equal for small values of r. Thus as del t -> 0, d is concave. But even
% for large values of r, the original function seems to be concave albeit
% it differs from the approximated one. Also, both attain their maximum at
% the same angle. 
