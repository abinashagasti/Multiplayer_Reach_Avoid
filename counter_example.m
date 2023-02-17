clear
clc

N=2;M=2;
a=[100 2;2 -1000];
f=a';f=a(:);
A2=[1 1 0 0;0 0 1 1];
A1=[1 0 1 0;0 1 0 1];
b=[1;1];

x_strict=linprog(-f,A1,b,A2,b,zeros(4,1));x_strict=reshape(x_strict,[N,M])'
x_relax=linprog(-f,[A1;A2],[b;b],[],[],zeros(4,1));x_relax=reshape(x_relax,[N,M])'