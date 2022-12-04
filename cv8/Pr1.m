clc;
clear all;
close all;

syms kk r r1 r2 p A B

sigmaR = A + B/r^2
sigmaT = A - B/r^2
sigmaA = (sigmaR + sigmaT)/2

rc1 = subs(sigmaR,r,r1) == -p
rc2 = subs(sigmaR,r,r2) == 0

r1 = 50;
r2 = 75;
p = 200;
kk = 1.5;

[Q Res] = equationsToMatrix([rc1 rc2], [A B])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
A = vpa(sol(1,1),3)
B = vpa(sol(2,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)
sigmaA = subs(sigmaA)

fplot(sigmaR,[r1, r2],'color','red')
hold on
fplot(sigmaT,[r1, r2],'color','blue')
fplot(sigmaA,[r1, r2],'color','black')

sigmaMax = subs(sigmaT,r,r1) - subs(sigmaR,r,r1)

Re = kk*sigmaMax


