clc;
clear all;
close all;


syms Ra Rb q x c bet k E Iy G S

T1 = q*x - Ra;
Mo1 = -q*x^2/2 + Ra*x;

T2 = -Ra+q*(c/2+x)-Rb;
Mo2 = Ra*(x+c/2) - q*(x+c/2)^2/2 + Rb*x;

dT1Ra = diff(T1,Ra);
dMo1Ra = diff(Mo1,Ra);
dT2Rb = diff(T2,Rb);
dMo2Rb = diff(Mo2,Rb);
dT1Rb = diff(T1,Rb);
dMo1Rb = diff(Mo1,Rb);
dT2Ra = diff(T2,Ra);
dMo2Ra = diff(Mo2,Ra);

wA = 0 ==1/(E*Iy)*(int(Mo1*dMo1Ra,x,0,c/2) + int(Mo2*dMo2Ra,x,0,c/2)) + bet/(G*S)*(int(T1*dT1Ra,x,0,c/2) + + int(T2*dT2Ra,x,0,c/2));
wB = -Rb/k == 1/(E*Iy)*(int(Mo1*dMo1Rb,x,0,c/2) + int(Mo2*dMo2Rb,x,0,c/2)) + bet/(G*S)*int(T2*dT2Rb,x,0,c/2);

[Q Res] = equationsToMatrix([wA, wB], [Ra, Rb]);

q = 10;
Re = 250;
kk = 2;
d2 = 20;
d1 = 30;
c = 500;
S = pi*(d1^2-d2^2)/4;
E = 1.43e5;
G = E/(2*(1+0.3));
bet = 1.89;
Iy = pi*(d1^4-d2^4)/64
k = 5000;

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res);
Ra = vpa(sol(1,1),3)
Rb = vpa(sol(2,1),3)







