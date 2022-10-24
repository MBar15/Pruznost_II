clear all
clc
close all

syms FA Iy E a x MA

Mo1 = FA*x + MA;
dMo1FA = diff(Mo1,FA);
dMo1MA = diff(Mo1,MA);

wA = 0.5 == 1/(E * Iy) * (int(Mo1*dMo1FA,x,0,a))
fiA = 0 == 1/(E * Iy) * (int(Mo1*dMo1MA,x,0,a))

[Q, Res] = equationsToMatrix([wA,fiA], [FA, MA])

d = 30;
a = 400;
E = 2.1e5;
Iy = (pi*d^4) / 64;
Re = 235;
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q,Res)
FA = vpa(sol(1,1),3);
MA = vpa(sol(2,1),3);

sigma = vpa(abs(MA*32 / (pi*d^3)),3);
kk = floor(Re/sigma);


display(FA)
display(MA)
%display(sigma)
%display(kk)