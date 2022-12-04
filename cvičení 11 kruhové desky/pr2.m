clc;
clear all;
close all;

syms F T D E h mu r r1 r2 Q1(r) Q11(r) w(r) C1 C2 

T = F/(2*pi*r)
D = E*h^3/(12*(1-mu^2))

ode1 = diff(Q1,r,1) == -T/D

iQ1(r) = dsolve(ode1)

ode2 = 1/r*diff(Q11,r,1) == iQ1

iQ2(r) = dsolve(ode2)
fi = iQ2/r
dfi = diff(fi,r)

sigmaR = E*h/(2*(1-mu^2))*(dfi+mu*fi/r)
sigmaT = E*h/(2*(1-mu^2))*(fi/r + mu*dfi)

% ode3 = diff(w,r,1) == fi
% w = dsolve(ode3)


E = 70e3;
mu = 0.3;
r1 = 100;
r2 = 150;
h = 10;
Re = 280;
kk = 1.5

% OP
%rce1 = subs(w,r,r2) == 0
rce2 = subs(sigmaR,r,r1) == 0
rce3 = subs(sigmaR,r,r2) == 0

[Q Res] = equationsToMatrix([rce2 rce3], [C1 C2])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
C1 = vpa(sol(1,1),3)
C2 = vpa(sol(2,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)

fplot(subs(sigmaR,F,1),[r1 r2],'color','red')
hold on
fplot(subs(sigmaT,F,1),[r1 r2],'color','blue')

sigmaMax = subs(sigmaT,r,r1)

rce4 = sigmaMax == Re/kk
F = vpa(solve(rce4),4)

