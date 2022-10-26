clc;
clear all;
close all;

syms kk Re D d F R q beta alfa fi x

q = 2;
beta = pi/4;
alfa = pi/3;
R = 500;
F = 800;
Re = 250;
D1 = 40;
d1 = 32;
Wo = pi/(32*D1) * (D1^4-d1^4);

qy1 = q*cos(beta)
qx1 = q*sin(beta)

Fy1 = F*sin(alfa)
Fx1 = F*cos(alfa)

a = R*sin(fi);
b = R - R*cos(fi);
c = R*cos(fi);
d = R + R*sin(fi);

Mo1 = -Fy1 * b - Fx1*a;
Mo2 = -Fy1*c - Fx1*d;
Mo3 = -Fx1*2*R + Fy1*x -qy1*x^2/2;



Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);

fplot(Mo1,[0,pi/2],'color','red')
hold on
fplot(Mo2,[0,pi/2],'color','blue')
%fplot(Mo3,[0,2*R],'color','black')


sigmax = qx1*2*R / (pi*(D1^2-d1^2)/4)
Momax = vpa(abs(subs(Mo3,x,2*R)),3)
sigmao = vpa(Momax / Wo,3)

sigmamax = sigmax +  sigmao;

kk = Re/sigmamax
