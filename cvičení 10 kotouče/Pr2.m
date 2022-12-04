clc;
clear all;
close all;

syms A1 B1 A2 B2 r r1 r2 r3 ps rho1 rho2 OM mu1 mu2 E1 E2 delt

% OM = omega^2
sigmaR1 = A1 + B1/r^2 - (3 + mu1)/8*rho1*OM*r^2
sigmaT1 = A1 - B1/r^2 - (1 + 3*mu1)/8*rho1*OM*r^2

sigmaR2 = A2 + B2/r^2 - (3 + mu2)/8*rho2*OM*r^2
sigmaT2 = A2 - B2/r^2 - (1 + 3*mu2)/8*rho2*OM*r^2

rc1 = subs(sigmaR1,r,r1) == 0
rc2 = subs(sigmaR1,r,r2) == -ps
rc3 = subs(sigmaR2,r,r2) == -ps
rc4 = subs(sigmaR2,r,r3) == 0

sigmaR1R2 = subs(sigmaR1,r,r2);
sigmaR2R2 = subs(sigmaR2,r,r2);
sigmaT1R2 = subs(sigmaT1,r,r2);
sigmaT2R2 = subs(sigmaT2,r,r2);

u12 = r2/E1*(sigmaT1R2 - mu1*(sigmaR1R2))
u22 = r2/E2*(sigmaT2R2 - mu2*(sigmaR2R2))

rc5 = delt == -u12 + u22

E1 = 210e9;
E2 = 100e9;
rho1 = 7850;
rho2 = 4540;
mu1 = 0.3;
mu2 = 0.3;
r1 = 0.1;
r2 = 0.16;
r3 = 0.25;
delt = 0.03e-3;

[Q Res] = equationsToMatrix([rc1 rc2 rc3 rc4 rc5], [A1 B1 A2 B2 ps])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)
A1 = vpa(sol(1,1),3)
B1 = vpa(sol(2,1),3)
A2 = vpa(sol(3,1),3)
B2 = vpa(sol(4,1),3)
ps = vpa(sol(5,1),3)

% uvolnovaci otacky jsou charakterizovany ze tlak vymizi
rc6 = -ps == 0
OM = solve(rc6,OM)
omega = OM^0.5
n = vpa(omega/(2*pi))


A1 = subs(A1)
A2 = subs(A2)
B1 = subs(B1)
B2 = subs(B2)
ps = subs(ps)

sigmaR1 = subs(sigmaR1)
sigmaR2 = subs(sigmaR2)
sigmaT1 = subs(sigmaT1)
sigmaT2 = subs(sigmaT2)

fplot(sigmaR1,[r1 r2],'color','red')
hold on
fplot(sigmaR2,[r2 r3],'color','black')
fplot(sigmaT1,[r1 r2],'color','green')
fplot(sigmaT2,[r2 r3],'color','blue')
legend('sigmaR1','sigmaR2','sigmaT1','sigmaT2')

