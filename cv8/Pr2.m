clc;
clear all;
close all;

syms kk h Re r r1 r2 A B p t

sigmaR = A + B/r^2
sigmaT = A - B/r^2
sigmaA = (sigmaR + sigmaT)/2

rc1 = subs(sigmaR,r,r1) == 0
rc2 = subs(sigmaR,r,r1+t) == -p

r1 = 4000;
h = 4000;
rho = 1000;
g = 9.81;
p = rho*h*g*1e-6
Re = 600;
kk = 1.5


[Q Res] = equationsToMatrix([rc1 rc2], [A B])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
A = vpa(sol(1,1),3)
B = vpa(sol(2,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)
sigmaA = subs(sigmaA)

fplot(subs(sigmaR,t,1000), [r1, r1+1000], 'color','red')
hold on
fplot(subs(sigmaT,t,1000), [r1, r1+1000], 'color','blue')
fplot(subs(sigmaA,t,1000), [r1, r1+1000], 'color','black')

sigmaMax = subs(sigmaR,r,r1) - subs(sigmaT,r,r1);
rc2 = sigmaMax == Re/kk
t = vpa(solve(rc2,t),3)