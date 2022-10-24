clear all
clc
close all

syms FA F c v Iy E a x Md

Mo1 = FA*x + Md
dMo1 = diff(Mo1,Md)

Mo2 = FA*(c+x) - F*x + Md
dMo2 = diff(Mo2, Md)

Iy = 1/12 * a^4

fiA = 1/(E*Iy) * (int(Mo1*dMo1,x,0,c) + int(Mo2*dMo2,x,0,c))

a = 20
E = 2.1e5
c = 500
F = 100
v = 0.5
FA = 27.1
Md = 0

%wA = vpa(subs(wA),3)
%FA = vpa(solve(wA,FA),3)
fiA = vpa(subs(fiA) * 180 / pi,3)