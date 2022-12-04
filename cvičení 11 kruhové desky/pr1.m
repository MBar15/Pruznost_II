clc;
clear all;
close all;

syms Re mu r2 r p h  Q1(r) Q11(r) T D C1 C2 C3 w(r) r2 E

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


% OP
rce1 = C2 == 0
rce2 = subs(sigmaR,r,r2) == 0
rce3 = subs(w,r,r2) == 0

% zadano
Re = 400;
mu = 0.3;
r2 = 350;
p = 1;
h = 20;
E = 2.1e5;

[Q Res] = equationsToMatrix([rce1 rce2 rce3], [C1 C2 C3])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
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

sigmaMax = subs(sigmaR,r,000000.1)
kk = Re/sigmaMax

hold off

fplot(w,[0 r2], 'color','green')
wmax = subs(w,r,0.0000001)
