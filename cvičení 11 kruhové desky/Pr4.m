clc;
clear all;
close all;

syms T D E h r p mu Q1(r) Q11(r) w(r) C1 C2 C3 r2

T = p*r/2
D = E*h^3/(12*(1-mu^2))

ode1 = diff(Q1,r,1) == -T/D

iQ1(r) = dsolve(ode1)

ode2 = 1/r*diff(Q11,r,1) == iQ1
iQ2(r) = dsolve(ode2)

fi = iQ2(r)/r
dfi = diff(fi,r)

sigmaR = E*h/(2*(1-mu^2))*(dfi+mu*fi/r)
sigmaT = E*h/(2*(1-mu^2))*(fi/r + mu*dfi)

ode3 = diff(w,r,1) == fi
w = dsolve(ode3)


% OP
rce1 = subs(w,r,r2) == 0
rce2 = subs(fi,r,r2) == 0
rce3 = C2 == 0
E = 210e3;
mu = 0.35;
r2 = 350;
p = 1;
h =20;
Rmtah = 220;
Rmtlak = 450;

[Q Res] = equationsToMatrix([rce1 rce2 rce3], [C1 C2 C3])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q,Res)
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

sigmaMax = sigmaZ - subs(sigmaR,r,r2)

kk = Rmtah/sigmaMax

hold off

fplot(w,[0 r2], 'color','green')
wmax = subs(w,r,0.0000001)
