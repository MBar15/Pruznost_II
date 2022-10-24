clc;
clear all;
close all;

syms q b h l F E Iy beta1 G S x

Mo1 = -F*x-q*x^2/2;
dMo1F = diff(Mo1,F);
T = F+q*x
dT = diff(T,F)
wA = 5 == 1/(E*Iy) * int(Mo1*dMo1F,0,l) + beta1/(G*S)*(int(T*dT,x,0,l));

b = 20;
h = 30;
l = 500;
F = 500;
E = 143e3 %2.1e5
Iy = 1/12*b*h^3;
beta1 = 1 %2;
G = 2500%E/(2*(1+0.3));
S = b*h;

wA = subs(wA)

q = vpa(solve(wA,q),3)










