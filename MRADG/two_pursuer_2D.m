clear
clc

syms xc1 xc2 yc1 yc2 rc1 rc2 x y

C1=(x-xc1)*(x-xc1)+(y-yc1)*(y-yc1)-rc1*rc1;
C2=(x-xc2)*(x-xc2)+(y-yc2)*(y-yc2)-rc2*rc2;

L=simplify(C1-C2,1000);
ytemp=simplify(solve(L,y));
xint=simplify(solve(subs(C1,y,ytemp),x),1000);
[N,D]=numden(xint-xc1/2-xc2/2);
yint=simplify(subs(ytemp,x,xint),1000);

%%

syms A B C
[xint1,A]=subexpr(xint,A);
[xint2,B]=subexpr(xint1,B);
[xsol,C]=subexpr(xint2,C);
[yint1,A]=subexpr(yint,A);
[yint2,B]=subexpr(yint1,B);
[ysol,C]=subexpr(yint2,C);

%%

% xsol1=xc1/2+xc2/2+((rc1^2-rc2^2)*(xc2-xc1))/(2*(xc1-xc2)^2+2*(yc1-yc2)^2)+((yc1-yc2)*sqrt(4*rc1^2*rc2^2-((xc1-xc2)^2+(yc1-yc2)^2-rc1^2-rc2^2)^2))/(2*(xc1-xc2)^2+2*(yc1-yc2)^2);
% xsol2=xc1/2+xc2/2+((rc1^2-rc2^2)*(xc2-xc1))/(2*(xc1-xc2)^2+2*(yc1-yc2)^2)-((yc1-yc2)*sqrt(4*rc1^2*rc2^2-((xc1-xc2)^2+(yc1-yc2)^2-rc1^2-rc2^2)^2))/(2*(xc1-xc2)^2+2*(yc1-yc2)^2);
% 
% %simplify(xint(1)-xsol1,1000)
% 
% ysol1=yc1/2+yc2/2+((rc1^2-rc2^2)*(yc2-yc1))/(2*(xc1-xc2)^2+2*(yc1-yc2)^2)+((xc2-xc1)*sqrt(4*rc1^2*rc2^2-((xc1-xc2)^2+(yc1-yc2)^2-rc1^2-rc2^2)^2))/(2*(xc1-xc2)^2+2*(yc1-yc2)^2);
% ysol2=yc1/2+yc2/2+((rc1^2-rc2^2)*(yc2-yc1))/(2*(xc1-xc2)^2+2*(yc1-yc2)^2)-((xc2-xc1)*sqrt(4*rc1^2*rc2^2-((xc1-xc2)^2+(yc1-yc2)^2-rc1^2-rc2^2)^2))/(2*(xc1-xc2)^2+2*(yc1-yc2)^2);
% 
% %simplify(ysol(1)-ysol1,1000)
% 
% xcn1=-2;ycn1=-1;xcn2=1;ycn2=2;
% rcn1=3;rcn2=2;
% 
% xsol1=(xc1+xc2)/2+(B*(xc2-xc1))/C+(A*(yc1-yc2))/(2*C);
% xsol2=(xc1+xc2)/2+(B*(xc2-xc1))/C-(A*(yc1-yc2))/(2*C);
% ysol1=(yc1+yc2)/2+(B*(yc2-yc1))/C-(A*(xc1-xc2))/(2*C);
% ysol2=(yc1+yc2)/2+(B*(yc2-yc1))/C+(A*(xc1-xc2))/(2*C);
% 
% Aval=double(subs(A,[xc1,xc2,yc1,yc2,rc1,rc2],[xcn1,xcn2,ycn1,ycn2,rcn1,rcn2]));
% Bval=double(subs(B,[xc1,xc2,yc1,yc2,rc1,rc2],[xcn1,xcn2,ycn1,ycn2,rcn1,rcn2]));
% Cval=double(subs(C,[xc1,xc2,yc1,yc2,rc1,rc2],[xcn1,xcn2,ycn1,ycn2,rcn1,rcn2]));
% 
% xsol1n=subs(xsol1,[xc1,xc2,yc1,yc2,rc1,rc2,A,B,C],[xcn1,xcn2,ycn1,ycn2,rcn1,rcn2,Aval,Bval,Cval]);
% ysoln=subs(ysol,[xc1,xc2,yc1,yc2,rc1,rc2],[xcn1,xcn2,ycn1,ycn2,rcn1,rcn2]);
% 
% xsoln(1)=subs(xsoln(1),A,Aval)


% subs(xsol2,[xc1,xc2,yc1,yc2,rc1,rc2],[xcn1,xcn2,ycn1,ycn2,rcn1,rcn2])
% subs(ysol2,[xc1,xc2,yc1,yc2,rc1,rc2],[xcn1,xcn2,ycn1,ycn2,rcn1,rcn2])

%%

Vguess=xsol.^2+ysol.^2;

syms Rc1 Rc2
Rc1=sqrt(xc1^2+yc1^2);
Rc2=sqrt(xc2^2+yc2^2);

Vcomp=0.25*(Rc1^2+Rc2^2)+0.5*(xc1*xc2+yc1*yc2)+B^2/C+0.25*A^2/C+(B/C)*(Rc2^2-Rc1^2)+(A/C)*(xc2*yc1-xc1*yc2);
Vgrad=gradient(Vcomp,[xc1,xc2,yc1,yc2,rc1,rc2]);

%%
clear
clc

syms xc1 xc2 yc1 yc2 rc1 rc2 A B C Rc1 Rc2 xe ye x1 y1 x2 y2 alp1 alp2


xterm1 = -(4*(xc1-xc2)*B^2)/C^2;
xterm2 = (4*(xc1-xc2)*(Rc1^2-Rc2^2)*B)/C^2;
xterm3 = xc1+xc2;
xterm4 = (4*(xc1-xc2)*(C-2*rc1*rc2)*(xc1*yc2-xc2*yc1))/(A*C);
xterm5 = (4*A*(xc1-xc2)*(xc1*yc2-xc2*yc1))/C^2;
xterm61 = -4*B*(xc1/C);
xterm62 = 4*B*(xc2/C);
xterm7 = -2*((xc1-xc2)*(C-2*rc1*rc2))/C;
xterm81 = -2*A*(yc2/C);
xterm82 = 2*A*(yc1/C);
xterm9 = -(xc1-xc2)*(A^2/C^2);

yterm1 = -(4*(yc1-yc2)*B^2)/C^2;
yterm2 = (4*(yc1-yc2)*(Rc1^2-Rc2^2)*B)/C^2;
yterm3 = yc1+yc2;
yterm4 = (4*(yc1-yc2)*(C-2*rc1*rc2)*(xc1*yc2-xc2*yc1))/(A*C);
yterm5 = (4*A*(yc1-xc2)*(xc1*yc2-xc2*yc1))/C^2;
yterm61 = -4*B*(yc1/C);
yterm62 = 4*B*(yc2/C);
yterm7 = -2*((yc1-yc2)*(C-2*rc1*rc2))/C;
yterm81 = 2*A*(xc2/C);
yterm82 = -2*A*(xc1/C);
yterm9 = -(yc1-yc2)*(A^2/C^2);

rterm11 = 4*rc2*(C-2*rc1*rc2);
rterm12 = 4*rc1*(C-2*rc1*rc2);
rterm21 = 8*(xc2*yc1-xc1*yc2)*((rc2*C-rc1*B)/A);
rterm22 = 8*(xc2*yc1-xc1*yc2)*((rc1*C+rc2*B)/A);
rterm31 = 8*rc1*(Rc2^2-Rc1^2);
rterm32 = -8*rc2*(Rc2^2-Rc1^2);
rterm41 = 8*rc1*rc2^2;
rterm42 = 8*rc2*rc1^2;
rterm51 = 12*rc1*B;
rterm52 = -12*rc2*B;

Vxc1 = 0.5*(xterm1+xterm2+xterm3+xterm4+xterm5+xterm61+ ...
    xterm7+xterm81+xterm9);
Vyc1 = 0.5*(yterm1+yterm2+yterm3+yterm4+yterm5+yterm61+ ...
    yterm7+yterm81+yterm9);
Vxc2 = 0.5*(-xterm1-xterm2+xterm3-xterm4-xterm5+xterm62 ...
    -xterm7+xterm82-xterm9);
Vyc2 = 0.5*(-yterm1-yterm2+yterm3-yterm4-yterm5+yterm62 ...
    -yterm7+yterm82-yterm9);
Vrc1 = (0.25/C)*(rterm11+rterm21+rterm31+rterm41+rterm51);
Vrc2 = (0.25/C)*(rterm12+rterm22+rterm32+rterm42+rterm52);

xc1xe = 1/(1-alp1^2);
xc2xe = 1/(1-alp2^2);
rc1xe = (alp1/(1-alp1^2))*((xe-x1)/((xe-x1)^2+(ye-y1)^2));
rc2xe = (alp2/(1-alp2^2))*((xe-x2)/((xe-x2)^2+(ye-y2)^2));

yc1ye = 1/(1-alp1^2);
yc2ye = 1/(1-alp2^2);
rc1ye = (alp1/(1-alp1^2))*((ye-y1)/((xe-x1)^2+(ye-y1)^2));
rc2ye = (alp2/(1-alp2^2))*((ye-y2)/((xe-x2)^2+(ye-y2)^2));

xc1x1 = -alp1^2/(1-alp1^2);
rc1x1 = (alp1/(1-alp1^2))*((x1-xe)/((xe-x1)^2+(ye-y1)^2));

yc1y1 = -alp1^2/(1-alp1^2);
rc1y1 = (alp1/(1-alp1^2))*((y1-ye)/((xe-x1)^2+(ye-y1)^2));

xc2x2 = -alp2^2/(1-alp2^2);
rc2x2 = (alp2/(1-alp2^2))*((x2-xe)/((xe-x2)^2+(ye-y2)^2));

yc2y2 = -alp2^2/(1-alp2^2);
rc2y2 = (alp2/(1-alp2^2))*((y2-ye)/((xe-x2)^2+(ye-y2)^2));

Vxe = Vxc1*xc1xe+Vxc2*xc2xe+Vrc1*rc1xe+Vrc2*rc2xe;
Vye = Vyc1*yc1ye+Vyc2*yc2ye+Vrc1*rc1ye+Vrc2*rc2ye;
Vx1 = Vxc1*xc1x1+Vrc1*rc1x1;
Vy1 = Vyc1*yc1y1+Vrc1*rc1y1;
Vx2 = Vxc2*xc2x2+Vrc2*rc2x2;
Vy2 = Vyc2*yc2y2+Vrc2*rc2y2;

rhoe = sqrt(Vxe^2+Vye^2);
rho1 = sqrt(Vx1^2+Vy1^2);
rho2 = sqrt(Vx2^2+Vy2^2);

me2=-rhoe+rho1/alp1+rho2/alp2;

syms alpha beta gamma delta epsilon zeta eta theta
[me21,alpha]=subexpr(me2,alpha);
[me22,beta]=subexpr(me21,beta);
[me23,gamma]=subexpr(me22,gamma);
[me24,delta]=subexpr(me23,delta);
[me25,epsilon]=subexpr(me24,epsilon);
[me26,zeta]=subexpr(me25,zeta);
[me27,eta]=subexpr(me26,eta);
[me28,theta]=subexpr(me27,theta);

% alpha=8*rc1^2*rc2 - 12*B*rc2 + 4*rc1*(C - 2*rc1*rc2) + 8*rc2*(Rc1^2 - Rc2^2) - ((8*xc1*yc2 - 8*xc2*yc1)*(B*rc2 + C*rc1))/A;
% beta=12*B*rc1 + 8*rc1*rc2^2 + 4*rc2*(C - 2*rc1*rc2) - 8*rc1*(Rc1^2 - Rc2^2) + ((8*xc1*yc2 - 8*xc2*yc1)*(B*rc1 - C*rc2))/A;
% gamma=xc1/2 + xc2/2 - (A^2*(xc1 - xc2))/(2*C^2) - ((xc1 - xc2)*(C - 2*rc1*rc2))/C - (B^2*(4*xc1 - 4*xc2))/(2*C^2) - (2*B*xc1)/C - (A*yc2)/C + (B*(4*xc1 - 4*xc2)*(Rc1^2 - Rc2^2))/(2*C^2) + (2*A*(xc1*yc2 - xc2*yc1)*(xc1 - xc2))/C^2 + ((xc1*yc2 - xc2*yc1)*(4*xc1 - 4*xc2)*(C - 2*rc1*rc2))/(2*A*C);
% delta=yc1/2 + yc2/2 - (A^2*(yc1 - yc2))/(2*C^2) - ((yc1 - yc2)*(C - 2*rc1*rc2))/C - (B^2*(4*yc1 - 4*yc2))/(2*C^2) + (A*xc2)/C - (2*B*yc1)/C + (B*(4*yc1 - 4*yc2)*(Rc1^2 - Rc2^2))/(2*C^2) - (2*A*(xc1*yc2 - xc2*yc1)*(xc2 - yc1))/C^2 + ((xc1*yc2 - xc2*yc1)*(4*yc1 - 4*yc2)*(C - 2*rc1*rc2))/(2*A*C);
% epsilon=yc1/2 + yc2/2 + (A^2*(yc1 - yc2))/(2*C^2) + ((yc1 - yc2)*(C - 2*rc1*rc2))/C + (B^2*(4*yc1 - 4*yc2))/(2*C^2) - (A*xc1)/C + (2*B*yc2)/C - (B*(4*yc1 - 4*yc2)*(Rc1^2 - Rc2^2))/(2*C^2) + (2*A*(xc1*yc2 - xc2*yc1)*(xc2 - yc1))/C^2 - ((xc1*yc2 - xc2*yc1)*(4*yc1 - 4*yc2)*(C - 2*rc1*rc2))/(2*A*C);
% zeta=xc1/2 + xc2/2 + (A^2*(xc1 - xc2))/(2*C^2) + ((xc1 - xc2)*(C - 2*rc1*rc2))/C + (B^2*(4*xc1 - 4*xc2))/(2*C^2) + (A*yc1)/C + (2*B*xc2)/C - (B*(4*xc1 - 4*xc2)*(Rc1^2 - Rc2^2))/(2*C^2) - (2*A*(xc1*yc2 - xc2*yc1)*(xc1 - xc2))/C^2 - ((xc1*yc2 - xc2*yc1)*(4*xc1 - 4*xc2)*(C - 2*rc1*rc2))/(2*A*C);
% eta=1/(alp1^2 - 1);
% theta=1/((x2 - xe)^2 + (y2 - ye)^2);

% alpha=4*C*Vrc2;
% beta=4*C*Vrc1;
% gamma=Vxc1;         
% delta=Vyc1;
% epsilon=Vyc2;
% zeta=Vxc2;

