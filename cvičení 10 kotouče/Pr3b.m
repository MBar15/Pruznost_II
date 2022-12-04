clc;
clear all;
close all;

syms r1 r2 r A1 B1 mu1 mu2 rho1 rho2 omega r3 OM A2 B2 ps E1 E2 delt

sigmaR1 = A1 + B1/r^2 - (3+mu1)/8*rho1*omega^2*r^2;
sigmaT1 = A1 - B1/r^2 - (1+3*mu1)/8*rho1*omega^2*r^2;

sigmaR2 = A2 + B2/r^2 - (3+mu2)/8*rho2*omega^2*r^2;
sigmaT2 = A2 - B2/r^2 - (1+3*mu2)/8*rho2*omega^2*r^2;

rce1 = subs(sigmaR1,r,r1) == 0;
rce2 = subs(sigmaR1,r,r2) == -ps;
rce3 = subs(sigmaR2,r,r2) == -ps;
rce4 = subs(sigmaR2,r,r3) == 0;

sigmaR1R2 = subs(sigmaR1,r,r2);
sigmaR2R2 = subs(sigmaR2,r,r2);

sigmaT1R2 = subs(sigmaT1,r,r2);
sigmaT2R2 = subs(sigmaT2,r,r2);

u12 = r2/E1*(sigmaT1R2 - mu1*sigmaR1R2)
u22 = r2/E2*(sigmaT2R2 - mu2*sigmaR2R2)

rce5 = delt == -u12 + u22

B1 = 0;
E1 = 100e9;
E2 = 70e9;
mu1 = 0.35
mu2 = 0.3
Re1 = 950e6;
Re2 = 280e6;
rho1 = 4540;
rho2 = 2800;
r2 = 0.2;
r3 = 0.3;
n = 5000/60;
omega = 2*pi*n
delt = 6.48e-4;

[Q Res] = equationsToMatrix([rce2 rce3 rce4 rce5], [A1 A2 B2 ps])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
A1 = vpa(sol(1,1),3)
A2 = vpa(sol(2,1),3)
B2 = vpa(sol(3,1),3)
ps = vpa(sol(4,1),3)

sigmaR1=subs(sigmaR1);
sigmaR2=subs(sigmaR2);
sigmaT1=subs(sigmaT1);
sigmaT2=subs(sigmaT2);

fplot(sigmaR1,[0 r2],'color','red')
hold on
fplot(sigmaR2,[r2 r3],'color','black')
fplot(sigmaT1,[0 r2],'color','blue')
fplot(sigmaT2,[r2 r3],'color','green')