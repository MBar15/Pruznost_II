clc;
clear all;
close all;

syms a Re d Fd Md x E Iy


% výpočet Fd = 783.0 a Md = -1.57e+5
Mo1 = Fd * x + Md;
dMo1Fd = diff(Mo1,Fd);
dMo1Md = diff(Mo1,Md);

wFd = 0.5 == 1/(E*Iy) * (int(Mo1*dMo1Fd,x,0,a));
fiMd = 0 ==  1/(E*Iy) * (int(Mo1*dMo1Md,x,0,a));

[Q Res] = equationsToMatrix([wFd, fiMd], [Fd, Md]);

a = 400;
d = 30;
Iy = (pi*d^4)/64;
E = 2.1e5;
Re = 235;

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res);

Fd = vpa(sol(1,1),3)
Md = vpa(sol(2,1),3)

% bezpečnost kk = 3,98 (3 zaokrouhleno dolů)
Momax = abs(Md);
Wo = (pi*d^3)/32;
sigma = Momax/Wo;
kk = vpa(Re/sigma,3)