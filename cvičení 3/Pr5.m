clc;
clear all;
close all;

syms d D l betaocel q c Ra Rb E Iy x G S k
 
% 2 Staticky neurčité

Mo1 = Ra*x - q*x^2/2  
T1 = q*x - Ra

Mo2 = Ra*(c/2+x) - q*(c/2+x)^2/2 + Rb * x
T2 = -Ra -Rb + q*(c/2+x)

dMo1Ra = diff(Mo1,Ra);
dT1Ra = diff(T1,Ra);
dMo2Ra = diff(Mo2,Ra);
dT2Ra = diff(T2,Ra);

dMo1Rb = diff(Mo1,Rb);
dT1Rb = diff(T1,Rb)
dMo2Rb = diff(Mo2,Rb);
dT2Rb = diff(T2,Rb)

wA = 0 == 1/(E*Iy) * (int(Mo1*dMo1Ra,x,0,c/2) + int(Mo2*dMo2Ra,x,0,c/2)) + betaocel/(G*S)*(int(T1*dT1Ra,x,0,c/2) + int(T2*dT2Ra,x,0,c/2));
wB = -Rb/k == 1/(E*Iy) * (int(Mo1*dMo1Rb,x,0,c/2) + int(Mo2*dMo2Rb,x,0,c/2)) + betaocel/(G*S)*(int(T1*dT1Rb,x,0,c/2) + int(T2*dT2Rb,x,0,c/2));

q = 10;
k = 5000;
d = 20;
D = 30;
c = 500
betaocel = 1.89;
Iy = pi*(D^4-d^4)/64;
E = 1.43e5;
G = E/(2*(1+0.3));
S = pi*(D^2-d^2)/4;


[Q Res] = equationsToMatrix([wA,wB], [Ra, Rb]);
Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res);
Ra = vpa(sol(1,1),3)
Rb = vpa(sol(2,1),3)













