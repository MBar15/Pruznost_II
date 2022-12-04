clc;
clear all;
close all;

syms A1 B1 A2 B2 r r1 r2 r3 delt p ps E mu %sigmaT1

sigmaR1 = A1 + B1/r^2
sigmaR2 = A2 + B2/r^2
sigmaT1 = A1 - B1/r^2
sigmaT2 = A2 - B2/r^2

rce1 = subs(sigmaR1,r,r1) == -p;
rce2 = subs(sigmaR1,r,r2) == -ps;
rce3 = subs(sigmaR2,r,r2) == -ps;
rce4 = subs(sigmaR2,r,r3) == 0;

sigmaR1R2 = subs(sigmaR1,r,r2)
sigmaR2R2 = subs(sigmaR2,r,r2)

sigmaT1R2 = A1 - B1/r2^2
sigmaT2R2 = A2 - B2/r2^2

sigmaZ1R2 = (sigmaT1R2 + sigmaR1R2)/2;
sigmaZ2R2 = (sigmaT2R2 + sigmaR2R2)/2;

u12 = r2/E*(sigmaT1R2 - mu*(sigmaR1R2 + sigmaZ1R2))
u22 = r2/E*(sigmaT2R2 - mu*(sigmaR2R2 + sigmaZ2R2))

rce5 = delt == -u12 + u22

r1 =200;
r2 = 240;
r3 = 320;
Re = 300;
E = 2.1e5
p = 50
delt = 0.02;
mu = 0.3;
kk = 1.8;

rce1 = subs(rce1)
rce2 = subs(rce2)
rce3 = subs(rce3)
rce4 = subs(rce4)
rce5 = subs(rce5)

[Q Res] = equationsToMatrix([rce1 rce2 rce3 rce4 rce5], [A1 B1 A2 B2 ps])
sol = linsolve(Q, Res)

A1 = vpa(sol(1,1),3)
B1 = vpa(sol(2,1),3)
A2 = vpa(sol(3,1),3)
B2 = vpa(sol(4,1),3)
ps = vpa(sol(5,1),3)

sigmaR1 = subs(sigmaR1)
sigmaR2 = subs(sigmaR2)
sigmaT1 = subs(sigmaT1)
sigmaT2 = subs(sigmaT2)


fplot(sigmaR1, [r1 r2], 'color','red')
hold on
fplot(sigmaR2, [r2 r3], 'color','blue')
fplot(sigmaT1, [r1 r2], 'color','green')
fplot(sigmaT2, [r2 r3], 'color','black')

sigmared1 = subs(sigmaT1,r,r1) - subs(sigmaR1,r,r1)
sigmared2 = subs(sigmaT2,r,r2) - subs(sigmaR2,r,r2)
sigmared = max([sigmared1, sigmared2])

kk = Re/sigmared


