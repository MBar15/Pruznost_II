clc;
clear all;
close all;

syms Q11(r) Q21(r) T E mu h D r p C1 C2 C3 w(r) r2 X1 X2
 
% pocitam na spodni strane (hrane)
ode1 = diff(Q11,r,1) == -T/D
ode2 = diff(Q21,r,1) == -T/D


T = p*r/2
D = E*h^3/(12*(1-mu^2))
ode1 = subs(ode1,T,0)
ode2 = subs(ode2,T,F/(2*pi*r))

iQ11(r) = dsolve(ode1)
iQ11(r) = subs(iQ11(r),C1,X1)
iQ21(r) = dsolve(ode1)
iQ21(r) = subs(iQ12(r),C1,X2)

ode3 = 1/r*diff(Q12,r,1) == iQ11(r)
ode4 = 1/r*diff(Q22,r,1) == iQ12(r)
iQ21(r) = dsolve(ode2) 
iQ21(r) = dsolve(ode2) 

fi = iQ2(r)/r
dfi = diff(fi,r)

sigmaR = E*h/(2*(1-mu^2))*(dfi+mu*fi/r)
sigmaT = E*h/(2*(1-mu^2))*(fi/r + mu*dfi)

ode3 = diff(w,r,1) == fi
w = dsolve(ode3)

% O.P.
rce1 = C2 == 0
rce2 = subs(fi,r,r2) == 0
rce3 = subs(w,r,r2) == 0

E = 100e3;
r2 = 350;
mu = 0.35;
p = 1;
Re = 400;
h = 20;
RmTlak = 450;
RmTah = 220;

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
sigmaZ = 0
w = subs(w)

fplot(sigmaR,[0 r2], 'color','red')
hold on
fplot(sigmaT,[0 r2], 'color','green')
fplot(sigmaZ,[0 r2], 'color','blue')
legend('sigmaR','sigmaT','sigmaZ')

sigmared = -subs(sigmaR,r,r2)
kk = RmTah/sigmared

hold off

fplot(w,[0 r2], 'color','green')
wmax = subs(w,r,0.0000001)



