clc;
clear all;
close all;

syms Fcx Mc Fc q x l Iy E 

Fc = q*l/2

Mo1 = Fcx*x - Mc;
Mo2 = Fcx*2*l + Fc*x - q*x^2/2 - Mc;

dMo1Fcx = diff(Mo1,Fcx);
dMo2Fcx = diff(Mo2,Fcx);

dMo1Mc = diff(Mo1,Mc)
dMo2Mc = diff(Mo2,Mc)

fiC = 0 == 1/(E*Iy) * (int(Mo1*dMo1Mc,x,0,2*l) + int(Mo2*dMo2Mc,x,0,l/2));
uC = 0 == 1/(E*Iy) * (int(Mo1*dMo1Fcx,x,0,2*l) + int(Mo2*dMo2Fcx,x,0,l/2));

[Q, Res] = equationsToMatrix([fiC, uC], [Mc, Fcx]);

d = 20;
l = 500;
Re = 270;
q = 2;

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res);

Mc = vpa(sol(1,1),3)
Fcx = vpa(sol(2,1),3)

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);

fplot(Mo1,[0,2*l],'color','red');
hold on
fplot(Mo2,[0,l/2],'color','blue');

Momax = subs(Mo2,l/2)

sigma = Momax*32/(pi*d^3);

kk = vpa(Re/sigma,3)
% kk = 5,09
