clc;
clear all;
close all;

syms F a b h R t q Re rho Ra Rb x phi y z E G bet Iy S C1 C3 C2

% i = 3 - 3 - 1 -1 = 2 S.N

N1 = -Ra;
T1 = q*x + F;
Mo1 = q*x^2/2 + F*x;



N2 = -Ra*cos(phi) + F*sin(phi) + q*a/3*sin(phi)
T2 = Ra*sin(phi) + F*cos(phi) + q*a/3*cos(phi)
%Mo2 = Ra*R*(1-cos(phi)) + F*R*sin(phi) + q*a/3*(a/6+ R*sin(phi))
Mo2 = int(T2*R,phi) + C2;
rc2 = subs(Mo2,phi,0) == subs(Mo1,x,a/3);
C2 = solve(rc2,C2);
Mo2 = subs(Mo2);

N3 = Ra
T3 = Rb - q*a/3 - F;
%Mo3 = -q*a/3*(-a/6 + x)^2/2 - F*(-a/3 + x) + Rb * x + Ra*2*R;
Mo3 = int(T3,x) + C3
rc3 = subs(Mo3,x,0) == subs(Mo2,phi,pi);
C3 = solve(rc3,C3);
Mo3 = subs(Mo3);


F = 8000;
a = 300;
b = 40;
h = 120;
R = 80;
t = 10;
q = 30;
Re = 400;
Ri = R - h/2;
E = 210000;
G = 80e3; %E/(2*(1+0.3));
Iy = 1/12*b*h^3 - 1/12*(b-2*t)*(h-2*t)^3;
bet = 1.2;

S = b*h - (h-2*t)*(b-2*t)

Rn = vpa(S/(int(int(1/rho,Ri,Ri+h),-b/2,b/2) - int(int(1/rho,Ri+t,Ri+h-t),-b/2+t,b/2-t)),4)

e = vpa(R - Rn,4)

dN1Ra = diff(N1,Ra);
dN2Ra = diff(N2,Ra);
dN3Ra = diff(N3,Ra);

dN1Rb = diff(N1,Rb);
dN2Rb = diff(N2,Rb);
dN3Rb = diff(N3,Rb);

dT1Ra = diff(T1,Ra);
dT2Ra = diff(T2,Ra);
dT3Ra = diff(T3,Ra);

dT1Rb = diff(T1,Rb);
dT2Rb = diff(T2,Rb);
dT3Rb = diff(T3,Rb);

dMo1Ra = diff(Mo1,Ra);
dMo2Ra = diff(Mo2,Ra);
dMo3Ra = diff(Mo3,Ra);

dMo1Rb = diff(Mo1,Rb);
dMo2Rb = diff(Mo2,Rb);
dMo3Rb = diff(Mo3,Rb);

wA = 0 == 1/(E*S)*(int(N1*dN1Ra,x,0,a/3) + int(N2*dN2Ra*R,phi,0,pi) + int(N3*dN3Ra,x,0,a)) + ...
    1/(E*Iy)*(int(Mo1*dMo1Ra,x,0,a/3) + int(Mo2*dMo2Ra*R,phi,0,pi) + int(Mo3*dMo3Ra,x,0,a)) + ...
    bet/(G*S)*(int(T1*dT1Ra,x,0,a/3) + int(T2*dT2Ra*R,phi,0,pi) + int(T3*dT3Ra,x,0,a))

uB = 0 == 1/(E*S)*(int(N1*dN1Rb,x,0,a/3) + int(N2*dN2Rb*R,phi,0,pi) + int(N3*dN3Rb,x,0,a)) + ...
    1/(E*Iy)*(int(Mo1*dMo1Rb,x,0,a/3) + int(Mo2*dMo2Rb*R,phi,0,pi) + int(Mo3*dMo3Rb,x,0,a)) + ...
    bet/(G*S)*(int(T1*dT1Rb,x,0,a/3) + int(T2*dT2Rb*R,phi,0,pi) + int(T3*dT3Rb,x,0,a))


[Q Res] = equationsToMatrix([wA,uB], [Ra, Rb])
Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q,Res);
Ra = vpa(sol(1,1),3)
Rb = vpa(sol(2,1),3)



%sigmaO = Mo2*y/(S*e*(Rn-z))

% pak nasleduje vypočet max ohybu pokud je to v zakřivené oblasti počítá se
% to podle teorie silně zak. prutů, jinak je to ohyb klasický pak se určí
% průběh sigma v průřezu 

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);

N1 = subs(N1);
N2 = subs(N2);
N3 = subs(N3);

%fplot(Mo2,[0, pi],'color','black')
%fplot(N2/S,[0, pi],'color','blue')

fi2Max = diff(Mo2,phi) == 0
fi2max = solve(fi2Max,phi)

Mo2Max = subs(Mo2,phi,fi2max(1,1))

sigmaO1 = Mo1*(h/2)/(Iy) + N1/S;
sigmaO2 = Mo2Max*z/(S*e*(Rn-z)) %+ subs(N2,phi,fi2max(1,1))/S
sigmaO3 = Mo3*(h/2)/(Iy) + N3/S;


fplot(sigmaO1,[0, a/3],'color','red')
hold on
fplot(sigmaO3,[0, a],'color','green')
hold off

z1 = Rn-R-h/2
z2 = Rn-Ri

%fplot(sigmaO2,[ -82.279, 37.72],'color','blue')

sigmaO2Max = subs(sigmaO2,z,37.2)
sigmaO2Max1 = subs(sigmaO2,z,50)

kk = Re/sigmaO2Max
