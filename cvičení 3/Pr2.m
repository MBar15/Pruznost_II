clc;
clear all;
close all;

syms z1 z h b y t

U1 = int(int(z1,z1,[z,h/2]),y,-b/2,b/2);
U1max = subs(U1,z,h/2-t)

U2 = int(int(z1,z1,[z,h/2-t]),y,-t/2,t/2);
U2max = subs(U2,z,t-h/2)

U3 = int(int(z1,z1,[z,-h/2]),y,-b/2,b/2);


b = 80;
h = 160;
t = 20
Iy = 1.9e7;
l = 500;
F = 5000;


U1 = subs(U1)
U1max = subs(U1max)
U2 = subs(U2)
U2max = subs(U2max)
U3 = subs(U3)

Tau1 = (F*U1)/(Iy*b)
Tau2 = (F*(U1max+U2))/(Iy*t)
Tau3 = (F*(U2max+U3))/(Iy*b)

fplot(Tau1,z,[h/2-t,h/2],'color','red')
hold on
fplot(Tau2,z,[-h/2+t,h/2-t],'color','green')
fplot(Tau3,z,[-h/2,-h/2+t],'color','blue')

Taumax = vpa(subs(Tau2,z,0),3)
Momax = abs(-F*l)
sigma = Momax/(Iy)*h/2