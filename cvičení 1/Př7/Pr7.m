clear all
clc
close all

syms Ma Fa Na F x Iy E a 

Fa = F/2;
Na = F/2;

Mo1 = Ma - Fa*x;
Mo2 = Ma - Fa*a - Na*x;

dMo1Ma = diff(Mo1,Ma);
dMo2Ma = diff(Mo2,Ma);

fiA = 0 == 1/(E*Iy) * (int(Mo1*dMo1Ma,x,0,a) + int(Mo2*dMo2Ma,x,0,a));

D = 40;
d = 36;
a = 600;
F = 300;
Re = 2/3 * 600;
Iy = pi/(64)*(D^4 - d^4) 
E = 2.1e5;

fiA = subs(fiA);

Ma = solve(fiA,Ma)

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);

fplot(Mo1,[0,a],'color','red');
hold on
fplot(Mo2,[0,a],'color','black');

Momax = subs(Mo1,x,0);

sigma = Momax*32*D/(D^4-d^4);
kk = vpa(Re/sigma,3)

% kk = 3,06
