clear
clc
close all

x1=[0,1,2];x2=[1,3,4];
r1=3;r2=2;
beta=sqrt((x2(1)-x1(1))^2+(x2(3)-x1(3))^2);
R=[-((x2(2)-x1(2))*(x2(1)-x1(1)))/(norm(x2-x1)*beta),((x2(1)-x1(1))^2+(x2(3)-x1(3))^2)/(norm(x2-x1)*beta),...
    -((x2(2)-x1(2))*(x2(3)-x1(3)))/(norm(x2-x1)*beta); (x2(3)-x1(3))/beta, 0, (x1(1)-x2(1))/beta;...
    (x2(1)-x1(1))/norm(x2-x1), (x2(2)-x1(2))/norm(x2-x1), (x2(3)-x1(3))/norm(x2-x1)];
R*(x2-x1)'
k=(norm(x2)^2+norm(x1)^2+r1^2-r2^2-2*x1*x2')/(2*norm(x2-x1))
z=(r1^2-r2^2+(1-2*k)*norm(x2-x1)^2)/(2*norm(x2-x1))