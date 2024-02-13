clear
clc
close all

syms xe ye ze xp yp zp alpha

xc=(xe-alpha^2*xp)/(1-alpha^2);
yc=(ye-alpha^2*yp)/(1-alpha^2);
zc=(ze-alpha^2*zp)/(1-alpha^2);
rc=(alpha/(1-alpha^2))*norm([xe,ye,ze]-[xp,yp,zp]);
Rc=norm([xc,yc,zc]);
Re=norm([xe,ye,ze]);
Rp=norm([xp,yp,zp]);

gradVP=[(1/(1-alpha^2))*(xc/Rc-(alpha^2/(1-alpha^2))*((xe-xp)/rc));
    (1/(1-alpha^2))*(yc/Rc-(alpha^2/(1-alpha^2))*((ye-yp)/rc));
    (1/(1-alpha^2))*(zc/Rc-(alpha^2/(1-alpha^2))*((ze-zp)/rc));
    (alpha^2/(1-alpha^2))*(-xc/Rc+(1/(1-alpha^2))*((xe-xp)/rc));
    (alpha^2/(1-alpha^2))*(-yc/Rc+(1/(1-alpha^2))*((ye-yp)/rc));
    (alpha^2/(1-alpha^2))*(-zc/Rc+(1/(1-alpha^2))*((ze-zp)/rc))];

gradVE=[(-1/alpha)*(xe/Re);
    (-1/alpha)*(ye/Re);
    (-1/alpha)*(ze/Re);
    xp/Rp;
    yp/Rp;
    zp/Rp];