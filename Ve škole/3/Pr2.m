clc;
clear all;
close all;

syms F b t h x z z1 y Iy

U1 = int(int(z1,[z, h/2]),y,[-b/2, b/2])
U1max = subs(U1,z,h/2-t)

U2 = int(int(z1,[z, h/2-t]),y,[-t/2, t/2]);
U2max = subs(U2,z,-h/2+t);

U3 = int(int(z1,[z,-h/2+t]),y,[-b/2, b/2]);
U3max = subs(U3,z,-h/2);

Mo1 = -F*x;
Tau1 = F*U1/(Iy*b);
Tau2 = F*(U1max+U2)/(Iy*t);        % skládaní 
Tau3 = F*(U1max+U2max+U3)/(Iy*b);


F = 5000;
c = 300;
b = 80;
h = 160;
t = 20;;
Iy = 1.9e7;
Tau1 = subs(Tau1);
Tau2 = subs(Tau2);
Tau3 = subs(Tau3);


fplot(Tau1,z, [h/2-t, h/2],'color','red')
hold on
fplot(Tau2,z, [-h/2+t, h/2-t],'color','green')
fplot(Tau3,z, [-h/2, -h/2+t],'color','blue')




