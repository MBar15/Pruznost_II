clc;
clear all;
close all;

syms Nc Tc phi Ra x Rbx Mc E Iy x1 phi1 C1 C2 C3 C4 q R

Ra = 3/2*q*R;
Rbx = -1/2*q*R;


T1 = -Nc*sin(phi) + Tc*cos(phi);
T2 = -Nc*cos(phi) - Tc*sin(phi) - Ra*cos(phi);
T3 = Nc + Rbx + Ra;
T4 = Nc + Rbx + Ra -q*x;

Mo1 = int(T1*R,phi,0,phi1) + C1
C1 = Mc;
Mo1 = subs(Mo1)

Mo2 = int(T2*R,phi,0,phi1) + C2
rce2 = subs(Mo1,phi1,pi/2) == subs(Mo2,phi1,0);
C2 = solve(rce2,C2);
Mo2 = subs(Mo2);


Mo3 = int(T3,x,0,x1) + C3
rce3 = subs(Mo2,phi1,pi/2) == subs(Mo3,x1,0);
C3 = solve(rce3,C3);
Mo3 = subs(Mo3);

Mo4 = int(T4,x,0,x1) + C4
rce4 = subs(Mo3,x1,R) == subs(Mo4,x1,0);
C4 = solve(rce4,C4);
Mo4 = subs(Mo4);

dMo1Nc = diff(Mo1,Nc);
dMo1Tc = diff(Mo1,Tc);
dMo1Mc = diff(Mo1,Mc);

dMo2Nc = diff(Mo2,Nc);
dMo2Tc = diff(Mo2,Tc);
dMo2Mc = diff(Mo2,Mc);

dMo3Nc = diff(Mo3,Nc);
dMo3Tc = diff(Mo3,Tc);
dMo3Mc = diff(Mo3,Mc);

dMo4Nc = diff(Mo4,Nc);
dMo4Tc = diff(Mo4,Tc);
dMo4Mc = diff(Mo4,Mc);


wN = 0 == 1/(E*Iy) * (int(Mo1*dMo1Nc*R,phi1,0,pi/2) + int(Mo2*dMo2Nc*R,phi1,0,pi/2) + int(Mo3*dMo3Nc,x1,0,R) + int(Mo4*dMo4Nc,x1,0,R))
wT = 0 == 1/(E*Iy) * (int(Mo1*dMo1Tc*R,phi1,0,pi/2) + int(Mo2*dMo2Tc*R,phi1,0,pi/2) + int(Mo3*dMo3Tc,x1,0,R) + int(Mo4*dMo4Tc,x1,0,R))
fiC = 0 == 1/(E*Iy) * (int(Mo1*dMo1Mc*R,phi1,0,pi/2) + int(Mo2*dMo2Mc*R,phi1,0,pi/2) + int(Mo3*dMo3Mc,x1,0,R) + int(Mo4*dMo4Mc,x1,0,R))

[Q Res] = equationsToMatrix([wN, wT, fiC], [Nc, Tc, Mc]);

R = 800;
%Iy = pi*(d^4)/64;
q = 15;
kk = 1.8;
Re = 300;

Q = subs(Q);
Res = subs(Res);

sol = linsolve(Q,Res);
Nc = vpa(sol(1,1),3)
Tc = vpa(sol(2,1),3)
Mc = vpa(sol(3,1),3)

%{
fplot(Mo1,[0,pi/2],'color','red')
hold on
fplot(Mo2,[0,pi/2],'color','blue')
fplot(Mo3,[0,R],'color','black')
fplot(Mo4,[0,R],'color','green')

Momax = vpa(subs(Mo3m,x1,0),4);
N3 = -Tc;
S = pi*d^2/4
sigmaNorm = abs(Mo3*32/(pi*d^3));
d = vpa(solve(sigmaNorm,d),4)
N = F*cos(alfa);
%}

