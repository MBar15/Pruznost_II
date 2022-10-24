clear all
clc
close all

syms F MA d1 d2 E Iy a F x

Mo1 = (F/2)*x - MA;
dMo1 = diff(Mo1,MA);

Mo2 = (F/2)*x - (F/2)*a + MA;
dMo2 = diff(Mo2, MA);

fiA = 1/(E*Iy) * (int(Mo1*dMo1,x,0,a) + int(Mo2*dMo2,x,0,a));

d1 = 40;
d2 = 36;
a = 600;
E = 2.1e5;
Re = 400;
Iy = pi/(32*d1)*(d1^3-d2^3);
F = 300;


fiA = subs(fiA);
MA = vpa(solve(fiA),3);
sigma = vpa(abs(MA*32*d1 / (pi*(d1^4-d2^4))),3);

kk = floor(Re/sigma);


display(MA)
display(sigma)
display(kk)