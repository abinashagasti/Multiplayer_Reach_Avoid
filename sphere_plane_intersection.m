clear all;
syms alp bet gam x y z l1 l2 d rc a b c xc yc zc

% alp = a*xc+b*yc+c*zc+d
% bet^2= a^2+b^2+c^2
% gam^2= xc^2+yc^2+zc^2
%
l1=-1+(bet*bet*l2+2*(alp-d))/(2*alp);
fl2=4*rc^2*(l1+1)^2-(bet^2)*l2^2-4*(alp-d)*l2-4*gam^2;
sl2=simplify(solve(fl2,l2));
vfun=((-a*l2+2*l1*xc)^2+(-b*l2+2*l1*yc)^2+(-c*l2+2*l1*zc)^2)/((2*(l1+1))^2);

sf1=subs(vfun,l2,sl2(1));
sf2=subs(vfun,l2,sl2(2));

[sn1,d1]=numden(sf1);
[sn2,d2]=numden(sf2);
 

XC1=0.5*diff(diff(sn1,xc),xc);
YC1=0.5*diff(diff(sn1,yc),yc);
ZC1=0.5*diff(diff(sn1,zc),zc);
rsn1=simplify(sn1-XC1*xc^2-YC1*yc^2-ZC1*zc^2);
xC1=diff(rsn1,xc);
yC1=diff(rsn1,yc);
zC1=diff(rsn1,zc);
rsn11=simplify(rsn1-xC1*xc-yC1*yc-zC1*zc);
revfun1=(1/d1)*((rsn11*bet^2/(a^2+b^2+c^2))+gam^2*XC1+xC1*(alp-d)/a);

XC2=0.5*diff(diff(sn2,xc),xc);
YC2=0.5*diff(diff(sn2,yc),yc);
ZC2=0.5*diff(diff(sn2,zc),zc);
rsn2=simplify(sn2-XC2*xc^2-YC2*yc^2-ZC2*zc^2);
xC2=diff(rsn2,xc);
yC2=diff(rsn2,yc);
zC2=diff(rsn2,zc);
rsn21=simplify(rsn2-xC2*xc-yC2*yc-zC2*zc);
revfun2=(1/d2)*((rsn21*bet^2/(a^2+b^2+c^2))+gam^2*XC2+xC2*(alp-d)/a);


