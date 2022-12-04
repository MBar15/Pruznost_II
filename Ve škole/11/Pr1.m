clc;
clear all;
close all;

syms Q1(r) Q11(r) T E mu h D r p C1 C2 C3 w(r) r2
 
% pocitam na spodni strane (hrane)
ode = diff(Q1,r,1) == -T/D

T = p*r/2
D = E*h^3/(12*(1-mu^2))
ode1 = subs(ode)

iQ1(r) = dsolve(ode1)


ode2 = 1/r*diff(Q11,r,1) == iQ1(r)
iQ2(r) = dsolve(ode2) 

fi = iQ2(r)/r
dfi = diff(fi,r)

sigmaR = E*h/(2*(1-mu^2))*(dfi+mu*fi/r)
sigmaT = E*h/(2*(1-mu^2))*(fi/r + mu*dfi)

ode3 = diff(w,r,1) == fi
w = dsolve(ode3)

% O.P.
rce1 = C2 == 0
rce2 = subs(sigmaR,r,r2) == 0
rce3 = subs(w,r,r2) == 0

E = 2.1e5;
r2 = 350;
mu = 0.3
p = 1;
Re = 400;
h = 20;

rce1 = subs(rce1)
rce2 = subs(rce2)
rce3 = subs(rce3)

[Q Res] = equationsToMatrix([rce1 rce2 rce3], [C1 C2 C3])
sol = linsolve(Q,Res)
Q = subs(Q)
Res = subs(Res)
C1 = vpa(sol(1,1),3)
C2 = vpa(sol(2,1),3)
C3 = vpa(sol(3,1),3)

sigmaR = subs(sigmaR)
sigmaT = subs(sigmaT)
w = subs(w)

fplot(sigmaR,[0 r2], 'color','red')
hold on
fplot(sigmaT,[0 r2], 'color','green')
legend('sigmaR','sigmaT')

sigmared = subs(sigmaT,r,0.00001)
kk = Re/sigmared

hold off

fplot(w,[0 r2], 'color','green')
wmax = subs(w,r,0.0000001)

