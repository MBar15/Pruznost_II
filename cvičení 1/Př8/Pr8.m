clc;
clear all;
close all;

syms F Mc Nc x a Ra E Iy d


Mo1 = -F/2*x + Mc;
Mo2 = -F/2*a/2 + Mc - Nc*x;
Mo3 = -F/2*(a/2 - x) + Mc - Nc*a - Ra*x;

dMo1Mc = diff(Mo1,Mc);
dMo2Mc = diff(Mo2,Mc);
dMo3Mc = diff(Mo3,Mc);

dMo1Nc = diff(Mo1,Nc);
dMo2Nc = diff(Mo2,Nc);
dMo3Nc = diff(Mo3,Nc);

fiC = 0 == 1/(E*Iy) * (int(Mo1*dMo1Mc,x,0,a/2) + int(Mo2*dMo2Mc,x,0,a) + int(Mo3*dMo3Mc,x,0,a/2));
uC = 0 == 1/(E*Iy) * (int(Mo1*dMo1Nc,x,0,a/2) + int(Mo2*dMo2Nc,x,0,a) + int(Mo3*dMo3Nc,x,0,a/2));

[Q Res] = equationsToMatrix([fiC, uC], [Mc, Nc]);

a = 600;
F = 500;
Re = 2/3*370;
kk = 1.8;
Ra = F/2;

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res);
Mc = vpa(sol(1,1),3)
Nc = vpa(sol(2,1),3)

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);

fplot(Mo1,[0,a/2],'color','red');
hold on
fplot(Mo2,[0,a],'color','blue');
fplot(Mo3,[0,a/2],'color','green');

Momax = subs(Mo1,x,0);

sigma = Momax*32/(pi*d^3) == Re/kk;

d = solve(sigma,d)

% d = 15,64 mm

