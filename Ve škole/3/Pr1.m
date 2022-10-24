clc;
clear all;
close all;


syms F x



Mo1 = -F*x;

F = 5000;
c = 300;
b = 80;
h = 160;

Momax = subs(Mo1,x,c);
Momax = abs(subs(Momax));
sigma = vpa(Momax*6/(b*h^2),3)

tau = 3/2*F/(b*h)




