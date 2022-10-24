clc;
clear all;
close all;

%{
i = 3 - (1+1+3) = 2 S. NEURČITÉ

3 ČÁSTI NOSNÍKU (3 Mo) 3 průběhy

NAPĚTÍ OD POSOUVAJÍCÍH SIL ZANEDBÁVÁM
%}

syms F Mk q x alfa E Iy Ra Rb a

% x <0,a> 
Mo1 = F*sind(alfa)*x - q*x^2/2;
Mk1 = 0;

% x <0,a> 
Mo2 = F*sind(alfa)*(a+x) - q*a*(x+a/2) + Ra*x;
Mk2 = -Mk;

% x <0,a> 
Mo3 = F*sind(alfa)*(a*2+x) - q*a*(a/2+a+x) + Ra*(a+x) + Rb*x;
Mk3 = -Mk - Mk;

dMo1Ra = diff(Mo1,Ra);
dMo2Ra = diff(Mo2,Ra);
dMo3Ra = diff(Mo3,Ra);


dMo1Rb = diff(Mo1,Rb);
dMo2Rb = diff(Mo2,Rb);
dMo3Rb = diff(Mo3,Rb);


wA = 0 == 1/(E*Iy) * (int(Mo1*dMo1Ra,x,0,a) + int(Mo2*dMo2Ra,x,0,a) + int(Mo3*dMo3Ra,x,0,a));
wB = 0 == 1/(E*Iy) * (int(Mo1*dMo1Rb,x,0,a) + int(Mo2*dMo2Rb,x,0,a) + int(Mo3*dMo3Rb,x,0,a));

E = 2.1e5;
Mk = 100e3;
d = 25;
Re = 300;
a = 300
F = 50e3;
q = 4;
alfa = 1;
Iy = pi*d^4/64;



[Q, Res] = equationsToMatrix([wA,wB], [Ra, Rb]);
Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res);

Ra = vpa(sol(1,1),3)
Rb = vpa(sol(2,1),3)

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);

Mk1= subs(Mk1)
Mk2= subs(Mk2)
Mk3= subs(Mk3)

N = F*cosd(alfa);
S = pi*d^2/4;

sigmaNorm1 = (N/S + abs(Mo1*32/(pi*d^3)))
sigmaNorm2 = (N/S + abs(Mo2*32/(pi*d^3)))
sigmaNorm3 = (N/S + abs(Mo3*32/(pi*d^3)))

Tau1 = Mk1*16/(pi*d^3);
Tau2 = Mk2*16/(pi*d^3);
Tau3 = Mk3*16/(pi*d^3);

sigmaRed1 = sqrt(sigmaNorm1^2 + 4*Tau1^2);
sigmaRed2 = sqrt(sigmaNorm2^2 + 4*Tau2^2);
sigmaRed3 = sqrt(sigmaNorm3^2 + 4*Tau3^2);

fplot(sigmaRed1,[0,a],'color','red')
hold on
fplot(sigmaRed2,[0,a],'color','blue')
fplot(sigmaRed3,[0,a],'color','black')

% 144,652

sigmaMax = vpa(subs(sigmaRed3,x,0),3)

kk = vpa(Re/sigmaMax,3) % 1,71











