clc;
clear all;
close all;

syms kk h Re r r1 r2 A B p t

sigmaR = A + B/r^2
sigmaT = A - B/r^2
sigmaA = (sigmaR + sigmaT)/2

rc1 = subs(sigmaR,r,r1) == 0
rc2 = subs(sigmaR,r,r2) == -p

r1 = 4;
h = 4000;
rho = 1000;
g = 9.81;
p = rho*h*g
Re = 600e6;
kk = 1.5


[Q Res] = equationsToMatrix([rc1 rc2], [A B])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
A = (p*r2^2)/(r1^2 - r2^2)
B = -(p*r1^2*r2^2)/(r1^2 - r2^2)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)
sigmaA = subs(sigmaA)

fplot(subs(sigmaR,r2,1000), [r1, r1+0.5], 'color','red')
hold on
fplot(subs(sigmaT,r2,1000), [r1, r1+0.5], 'color','blue')
fplot(subs(sigmaA,r2,1000), [r1, r1+0.5], 'color','black')

sigmaMax = subs(sigmaR,r,r1) - subs(sigmaT,r,r1);
rc2 = sigmaMax == Re/kk
r2 = vpa(solve(rc2,r2),3)