clc;
clear all;
close all;

%{
    4 useky
    i = 3-2-1 = 0 S.N
    vnitrne 3 S.N.
%}

syms Tc Nc Mc E Iy x1 x phi C1 C2 C3 C4 q R phi1

Ra = 3/2*q*R;
Rb = 0;
Rbx = -1/2*q*R;

N1 = Nc*cos(phi) + Tc*sin(phi);
T1 = Tc*cos(phi) - Nc*sin(phi);
Mo1 = int(T1*R,phi,0,phi1) + C1;
C1 = Mc;
Mo1 = subs(Mo1);

N2 = Tc*cos(phi) - Nc*sin(phi) - Ra*sin(phi)
T2 = -Tc*sin(phi) - Nc*cos(phi) - Ra*cos(phi)
Mo2 = int(T2*R,phi,0,phi1) + C2;
rdce2 = subs(Mo1,phi1,pi/2) == subs(Mo2,phi1,0);
C2 = solve(rdce2,C2)
Mo2 = subs(Mo2)

N3 = -Tc;
T3 = Rbx + Ra + Nc;
Mo3 = int(T3,x,0,x1) + C3;
rdce3 = subs(Mo3,x1,0) == subs(Mo2,phi1,pi/2);
C3 = solve(rdce3,C3);
Mo3 = subs(Mo3);

N4 = -Tc
T4 = Rbx + Ra + Nc - q*x
Mo4 = int(T4,x,0,x1) + C4;
rdce4 = subs(Mo4,x1,0) == subs(Mo3,x1,R);
C4 = solve(rdce4,C4);
Mo4 = subs(Mo4);

dMo1Nc = diff(Mo1,Nc);
dMo2Nc = diff(Mo2,Nc);
dMo3Nc = diff(Mo3,Nc);
dMo4Nc = diff(Mo4,Nc);

dMo1Tc = diff(Mo1,Tc);
dMo2Tc = diff(Mo2,Tc);
dMo3Tc = diff(Mo3,Tc);
dMo4Tc = diff(Mo4,Tc);

dMo1Mc = diff(Mo1,Mc);
dMo2Mc = diff(Mo2,Mc);
dMo3Mc = diff(Mo3,Mc);
dMo4Mc = diff(Mo4,Mc);


wNc = 0  == 1/(E*Iy) * (int(Mo1*dMo1Nc*R,phi1,0,pi/2) + int(Mo2*dMo2Nc*R,phi1,0,pi/2) + int(Mo3*dMo3Nc,x1,0,R) + int(Mo4*dMo4Nc,x1,0,R));
wTc = 0  == 1/(E*Iy) * (int(Mo1*dMo1Tc*R,phi1,0,pi/2) + int(Mo2*dMo2Tc*R,phi1,0,pi/2) + int(Mo3*dMo3Tc,x1,0,R) + int(Mo4*dMo4Tc,x1,0,R));
fiMc = 0 == 1/(E*Iy) * (int(Mo1*dMo1Mc*R,phi1,0,pi/2) + int(Mo2*dMo2Mc*R,phi1,0,pi/2) + int(Mo3*dMo3Mc,x1,0,R) + int(Mo4*dMo4Mc,x1,0,R));

E = 2.1e5;
q = 15;
R = 800;


[Q Res] = equationsToMatrix([wNc,wTc,fiMc], [Nc,Tc,Mc]);

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q,Res);
Nc = vpa(sol(1,1),3)
Tc = vpa(sol(2,1),3)
Mc = vpa(sol(3,1),3)





