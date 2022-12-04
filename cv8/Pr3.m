clc;
clear all;
close all;

syms r r1 r2 p1 p2 kk Re A B p3

sigmaR = A + B/r^2
sigmaT = A - B/r^2
sigmaA = (sigmaR+ sigmaT)/2

rc1 = subs(sigmaR,r,r1) == -p1
rc2 = subs(sigmaR,r,r2) == -p3 %-p2

r1 = 600;
r2 = 800;
p2 = 30;
rho = 1000;
g = 9.81;
h = 3000;
p3 = rho*g*h*1e-6
kk = 1.5;
Re = 300;

[Q Res] = equationsToMatrix([rc1 rc2], [A B])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
A = vpa(sol(1,1),3)
B = vpa(sol(2,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)
sigmaA = subs(sigmaA)

fplot(subs(sigmaR,p1,40),[r1 r2],'color','red')
hold on
fplot(subs(sigmaT,p1,40),[r1 r2],'color','blue')
fplot(subs(sigmaA,p1,40),[r1 r2],'color','black')

sigmaMax1 = subs(sigmaT,r,r1) - subs(sigmaR,r,r1)
sigmaMax2 = subs(sigmaR,r,r1) - subs(sigmaT,r,r1)

rc3 = Re/kk == sigmaMax1
rc4 = Re/kk == sigmaMax2

p11 = vpa(solve(rc3,p1),3) % platí
p12 = solve(rc4,p1) % nelze vytvořit podlat zaporný

