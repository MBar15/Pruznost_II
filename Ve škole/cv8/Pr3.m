clc;
clear all;
close all;

syms A B p r1 r2 r g h rho t p2

sigmaR = A+B/r^2
sigmaT = A-B/r^2
sigmaA = (sigmaR+sigmaT)/2

rc1 = subs(sigmaR,r,r1) == -p
rc2 = subs(sigmaR,r,r2) == -h*rho*g

[Q Res] = equationsToMatrix([rc1, rc2], [A, B])

r1 = 0.6;
r2 = 0.8;
p2 = 30e6
kk = 1.5
rho = 1000;
h = 3000;
g = 9.81;
Re = 300e6

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res)
A = vpa(sol(1,1),3)
B = vpa(sol(2,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)
sigmaA = subs(sigmaA)


fplot(subs(sigmaR,p,20e6), [r1,r2],'color','red')
hold on
fplot(subs(sigmaT,p,20e6), [r1,r2],'color','blue')
fplot(subs(sigmaA,p,20e6), [r1,r2],'color','green')

sigmared= subs(sigmaT,r,r1)-subs(sigmaR,r,r1)

rc3 = sigmared == Re/kk
p = vpa(solve(rc3,p),3)