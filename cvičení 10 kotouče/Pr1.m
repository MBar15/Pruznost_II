clc;
clear all;
close all;

syms r r1 r2 mu A B rho omega

sigmaR = A + B/r^2 - (3 + mu)/8*rho*omega^2*r^2
sigmaT = A - B/r^2 - (1 + 3*mu)/8*rho*omega^2*r^2

rce1 = subs(sigmaR,r,r1) == 0
rce2 = subs(sigmaR,r,r2) == 0

r1 = 15/2*1e-3;
r2 = 120/2*1e-3;
h = 1.21e-3;
rho = 1190;
mu = 0.3
n = 10000/60;
E = 850e6;
Re = 60e6;
omega = 2*pi*n;

[Q Res] = equationsToMatrix([rce1 rce2], [A B])
Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res);
A = vpa(sol(1,1),3)
B = vpa(sol(2,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)

fplot(sigmaR,[r1,r2],'color','red')
hold on
fplot(sigmaT,[r1,r2],'color','blue')
hold off

sigmaMax = subs(sigmaT,r,r1) - subs(sigmaR,r,r1)

kk = vpa(Re/sigmaMax,3)

u = r/E*(sigmaT-mu*sigmaR)*1e3
fplot(u,[r1 r2], 'color','red')