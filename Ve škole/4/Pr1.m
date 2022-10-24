clc;
clear all;
close all;

syms alfa F q x Ra Rb E Iy a kk 

Mo1 = F*sind(alfa)*x - q*x^2/2;
Mo2 = F*sind(alfa)*(a+x) - q*a*(a/2+x) + Ra*x;
Mo3 = F*sind(alfa)*(2*a+x) - q*a*(x+3*a/2) + Ra*(a+x) + Rb*x;

dMo1Ra = diff(Mo1,Ra);
dMo1Rb = diff(Mo1,Rb);
dMo2Ra = diff(Mo2,Ra);
dMo2Rb = diff(Mo2,Rb);
dMo3Ra = diff(Mo3,Ra);
dMo3Rb = diff(Mo3,Rb);

wA = 0 == 1/(E*Iy)*(int(Mo1*dMo1Ra,x,0,a) + int(Mo2*dMo2Ra,x,0,a) +int(Mo3*dMo3Ra,x,0,a));
wB = 0 == 1/(E*Iy)*(int(Mo1*dMo1Rb,x,0,a) + int(Mo2*dMo2Rb,x,0,a) +int(Mo3*dMo3Rb,x,0,a));

[Q Res] = equationsToMatrix([wA, wB], [Ra, Rb]);

q = 4;
Re = 300;
d = 25;
a = 300;
E = 2.1e5;
Iy = pi*(d^4)/64;
F = 50e3;
alfa = 1;

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q,Res);

Ra = vpa(sol(1,1),3)
Rb = vpa(sol(2,1),3)

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);

%{
fplot(Mo1,[0,a],'color','red')
hold on
fplot(Mo2,[0,a],'color','blue')
fplot(Mo3,[0,a],'color','black')
%}

xMax = diff(Mo1,x) == 0;
xMax = vpa(solve(xMax,x),3)
Momax = vpa(subs(Mo1,x,xMax),3)

%sigma = Momax*32/(pi*d^3) + F/(pi*d^2/4) + 4*(Mk) == Re/kk;
%kk = solve(sigma,kk)

N = F*cosd(alfa);
S = pi*d^2/4;
sigmaNorm1 = (N/S + abs(Mo1)*32/(pi*d^3))
sigmaNorm2 = (N/S + abs(Mo2)*32/(pi*d^3))
sigmaNorm3 = (N/S + abs(Mo3)*32/(pi*d^3))

M1 = 100e3;

Tau2 = M1*16/(pi*d^3)
Tau3 = (M1*2)*16/(pi*d^3)

sigmared1 = sqrt(sigmaNorm1^2);
sigmared2 = sqrt(sigmaNorm2^2 + 4*Tau2^2);
sigmared3 = sqrt(sigmaNorm3^2 + 4*Tau3^2);

fplot(sigmared1,[0,a],'color','red')
hold on
fplot(sigmared2,[0,a],'color','blue')
fplot(sigmared3,[0,a],'color','black')

kk = vpa(Re/subs(sigmared3,x,0),3)
