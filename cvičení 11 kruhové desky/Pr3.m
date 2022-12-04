clc;
clear all;
close all;

syms T T1 T2 F r r1 r2 mu E h Q1(r) Q2(r) w(r) C1 C2 C3 X1 X2 X3 O1 O2 O3 fi fi1 fi2

% dva úseky
T1 = 0;
T2 = F/(2*pi*r)

D = E*h^3/(12*(1-mu^2))

ode = diff(Q1,r,1) == T/D

ode1 = subs(ode,T,T1)
ode2 = subs(ode,T,T2)

iQ11 = dsolve(ode1)
iQ11 = subs(iQ11,C1,X1)

iQ12 = dsolve(ode2)
iQ12 = subs(iQ12,C1,O1)

ode3 = 1/r*diff(Q2,r,1) == iQ11
ode4 = 1/r*diff(Q2,r,1) == iQ12

iQ21  = dsolve(ode3)
iQ21  = subs(iQ21,C1,X2)

iQ22  = dsolve(ode4)
iQ22  = subs(iQ22,C1,O2)

fi1 = X1*r/2+X2/r
dfi1 = diff(fi1,r)
fi2 = (O2 + (O1*r^2)/2 - (r^2*((3*F)/2 - (3*F*mu^2)/2 + 3*F*log(r)*(mu^2 - 1)))/(E*h^3*pi))/r
dfi2 = diff(fi2,r)

ode5 = diff(w,r,1) == fi1
w1 = dsolve(ode5)
w1 = subs(w1,C1,X3)

ode6 = diff(w,r,1) == fi2
w2 = dsolve(ode5)
w2 = subs(w2,C1,O3)


% rovnice napětí
sigmaR1 = E*h/(2*(1-mu^2))*(dfi1+mu*fi1/r)
sigmaT1 = E*h/(2*(1-mu^2))*(fi1/r + mu*dfi1)

sigmaR2 = E*h/(2*(1-mu^2))*(dfi2+mu*fi2/r)
sigmaT2 = E*h/(2*(1-mu^2))*(fi2/r + mu*dfi2)


% OP
rc1 = subs(sigmaR2,r,r2) == 0
rc2 = subs(sigmaR1,r,r1) == subs(sigmaR2,r,r1)
rc3 = subs(w1,r,r1) == subs(w2,r,r1)
rc4 = subs(fi1,r,r1) == subs(fi2,r,r1)
rc5 = subs(w2,r,r2) == 0
X2 = 0

E = 70e3;
mu = 0.3;
r1 = 100;
r2 = 200;
h = 15;
Re = 280;
F = 28e3;

[Q Res] = equationsToMatrix([rc1 rc2 rc3 rc4 rc5 ], [X1 X3 O1 O2 O3])
Q= subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
X1 = vpa(sol(1,1),3)
X3 = vpa(sol(2,1),3)
O1 = vpa(sol(3,1),3)
O2 = vpa(sol(4,1),3)
O3 = vpa(sol(5,1),3)

sigmaR1 = subs(sigmaR1)
sigmaT1 = subs(sigmaT1)
sigmaR2 = subs(sigmaR2)
sigmaT2 = subs(sigmaT2)
fi1 = subs(fi1)
fi2 = subs(fi2)
dfi1 = subs(dfi1)
dfi2 = subs(dfi2)
w1 = subs(w1)
w2 = subs(w2)

fplot(w1,[0,r1],'color','blue')
hold on
fplot(w2,[r1,r2],'color','red')
hold off

wmax = subs(w1,r,0)

fplot(sigmaR1,[0,r1],'color','red')
hold on
fplot(sigmaT1,[0,r1],'color','blue')
fplot(sigmaR2,[r1,r2],'color','green')
fplot(sigmaT2,[r1,r2],'color','black')


sigmaMax = -subs(sigmaT1,r,r1)

kk = Re/sigmaMax
