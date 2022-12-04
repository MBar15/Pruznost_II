clc;
clear all;
close all;

syms A1 B1 A2 B2 r r1 r2 r3 delt p ps E mu kk Re 

sigmaR1 = A1;
sigmaT1 = A1;

sigmaR2 = A2 + B2/r^2
sigmaT2 = A2 - B2/r^2

rce1 = sigmaR1 == sigmaT1
rce2 = subs(sigmaR1,r,r2) == -ps
rce3 = subs(sigmaR2,r,r2) == -ps
rce4 = subs(sigmaR2,r,r3) == 0

sigmaR1R2 = subs(sigmaR1,r,r2)
sigmaR2R2 = subs(sigmaR2,r,r2)

sigmaT1R2 = A1;
sigmaT2R2 = subs(sigmaT2,r,r2)

sigmaZ2R2 = (sigmaT2R2 + sigmaR2R2)/2

u12 = r2/E*(sigmaT1R2 - mu*(sigmaR1R2))
u22 = r2/E*(sigmaT2R2 - mu*(sigmaR2R2 + sigmaZ2R2))

rce5 = delt == -u12 + u22
rce6 = subs(sigmaT2,r,r2) - subs(sigmaR2,r,r2) == Re/kk

r2 = 25;
r3 =  50;
Re = 250;
kk = 1.5;
mu = 0.3;
f = 0.1;
b = 50;
E = 2.1e5

rce1 = subs(rce1)
rce2 = subs(rce2)
rce3 = subs(rce3)
rce4 = subs(rce4)
rce5 = subs(rce5)
rce6 = subs(rce6)

[Q Res] = equationsToMatrix([rce2 rce3 rce4 rce5 rce6], [A1 A2 B2 ps delt])
sol = linsolve(Q, Res)

A1 = vpa(sol(1,1),3)
A2 = vpa(sol(2,1),3)
B2 = vpa(sol(3,1),3)
ps = vpa(sol(4,1),3)
delt = vpa(sol(5,1),3)

sigmaR1 = subs(sigmaR1)
sigmaR2 = subs(sigmaR2)
sigmaT1 = subs(sigmaT1)
sigmaT2 = subs(sigmaT2)
sigmaZ2R2 = subs(sigmaZ2R2)

fplot(sigmaR1, [0 r2], 'color','red')
hold on
fplot(sigmaR2, [r2 r3], 'color','blue')
fplot(sigmaT1, [0 r2], 'color','green')
fplot(sigmaT2, [r2 r3], 'color','black')

Mk = vpa(2*pi*r2^2*ps*f*b,3)

