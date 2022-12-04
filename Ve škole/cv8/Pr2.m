clc;
clear all;
close all;

syms A B p r1 r2 r g h rho t

sigmaR = A+B/r^2
sigmaT = A-B/r^2
sigmaA = (sigmaR+sigmaT)/2

rc1 = subs(sigmaR,r,r1) == 0
rc2 = subs(sigmaR,r,r1+t) == -h*rho*g

[Q Res] = equationsToMatrix([rc1, rc2], [A, B])

r1 = 4;
kk = 1.5
rho = 1000;
h = 4000;
g = 9.81;
Re = 600e6

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res)
A = vpa(sol(1,1),3)
B = vpa(sol(2,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)
sigmaA = subs(sigmaA)


fplot(subs(sigmaR,t,0.5), [r1,r1+0.5],'color','red')
hold on
fplot(subs(sigmaT,t,0.5),[r1,r1+0.5],'color','blue')
fplot(subs(sigmaA,t,0.5), [r1,r1+0.5],'color','green')

sigmared= subs(sigmaR,r,r1)-subs(sigmaT,r,r1)

rc3 = sigmared == Re/kk
t = vpa(solve(rc3,t),3)