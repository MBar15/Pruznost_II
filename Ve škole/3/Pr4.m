clc;
clear all;
close all;

syms F x q G E b h bet c Iy S


Mo1 = -F*x-q*x^2/2
T1 = F + q*x;
%Mo = -int(T1,x)
dMo1 = diff(Mo1,F);
dT1 = diff(T1,F);

wA = 5 == 1/(E*Iy)*int(Mo1*dMo1,x,0,c) + bet/(G*S)*int(T1*dT1,x,0,c);

E = 1.43e5 %2.1e5;
G = 2500  %E/(2*(1+0.3))
h = 30;
b = 20;
Iy = 1/12 * b*h^3;
bet = 1;
S = b*h;
c = 500;
F = 500;

wA = subs(wA);
q = vpa(solve(wA,q),3)




