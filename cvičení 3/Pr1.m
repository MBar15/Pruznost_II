clc;
clear all;
close all;

syms F x Iy U b h y z

Mo1 = -F*x;
T1 = F

Tau = (T1*U)*(Iy*b)

U = int(int(z,z,0,h/2),y,-b/2,b/2)

b = 80;
h = 160;
l = 300
Iy = 1/12*b*h^3;
F = 5000;

Mo1 = subs(Mo1);
T1 = subs(T1)

%fplot(Mo1,[0,l],'color','red');
%hold on

Momax = subs(Mo1,x,l)

U = subs(U);
Wo = 1/6*b*h^2;

sigma = vpa(abs(Momax) / Wo,3)
Tau = vpa((T1*U)/(Iy*b),3)
