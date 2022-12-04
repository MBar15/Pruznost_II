clc;
clear all;
close all;

syms omega rho A1 B1 r mu r1 r2

sigmaR1 = A1 + B1/r^2 - (3+mu)/8*rho*omega^2*r^2;
sigmaT1 = A1-B1/r^2 - (1+3*mu)/8*rho*omega^2*r^2;

rc1 = 0 == subs(sigmaR1,r,r1)
rc2 = 0 == subs(sigmaR1,r,r2)

r1 = 0.0075;
r2 = 0.06;
Re = 60e6;
E = 850e6;
rho = 1190;
mu = 0.3;
n = 10000;
omega = 2*pi*n/60;

[Q Res] = equationsToMatrix([rc1, rc2],[A1 B1])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
A1 = vpa(sol(1,1),3)
B1 = vpa(sol(2,1),3)

sigmaR1 = subs(sigmaR1)
sigmaT1 = subs(sigmaT1)

fplot(sigmaR1,[r1, r2],'color','red')
hold on
fplot(sigmaT1,[r1, r2],'color','blue')
hold off

sigmaRed = subs(sigmaT1,r1) - subs(sigmaR1,r1) 
kk = vpa(Re/sigmaRed)

u = r/E*(sigmaT1-mu*sigmaR1)*1e3
fplot(u,[r1 r2],'color','red')