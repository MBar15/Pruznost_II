clc;
clear all;
close all;

syms Q1(r) Q11(r) T E mu h D r p C1 C2 C3 w(r) r1 r2 F
 
% pocitam na spodni strane (hrane)
ode = diff(Q1,r,1) == -T/D

T = F/(2*pi*r)
D = E*h^3/(12*(1-mu^2))
ode1 = subs(ode)

iQ1(r) = dsolve(ode1)


ode2 = 1/r*diff(Q11,r,1) == iQ1(r)
iQ2(r) = dsolve(ode2) 

fi = iQ2(r)/r
dfi = diff(fi,r)

sigmaR = E*h/(2*(1-mu^2))*(dfi+mu*fi/r)
sigmaT = E*h/(2*(1-mu^2))*(fi/r + mu*dfi)

% O.P.
rce1 = subs(sigmaR,r,r1) == 0
rce2 = subs(sigmaR,r,r2) == 0

E =70e3;
r1 = 100
r2 = 150;
mu = 0.3
Re = 280;
h = 10;
kk = 1.5

rce1 = subs(rce1)
rce2 = subs(rce2)

[Q Res] = equationsToMatrix([rce1 rce2], [C1 C2])
sol = linsolve(Q,Res)
Q = subs(Q)
Res = subs(Res)
C1 = vpa(sol(1,1),3)
C2 = vpa(sol(2,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)


fplot(subs(sigmaR,F,1000),[r1 r2], 'color','red')
hold on
fplot(subs(sigmaT,F,1000),[r1 r2], 'color','green')
legend('sigmaR','sigmaT')

sigmared = subs(sigmaT,r,r1)
 
rce4 = kk == Re/sigmared
F = vpa(solve(rce4),4)










