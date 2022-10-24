clc;
clear all;
close all;

syms kk Re D d F R q beta alfa fi x Fp

q = 2;
beta = pi/4;
alfa = pi/3;
R = 500;
F = 800;
Re = 250;
D1 = 40;
d1 = 32;
Wo = pi/D1 * (D1^4-d1^4);

qy1 = q*cos(beta);
qx1 = q*sin(beta);

Fx1 = F*sin(alfa);
Fy1 = F*cos(alfa);

a = R*sin(fi);
b = R - R*cos(fi);
c = R*cos(fi);
d = R + R*sin(fi);

Mo1 = -Fx1 * b - Fy1*a;
Mo2 = -Fy1*c - Fx1*d;
Mo3 = -Fx1*2*R + Fy1*x -qy1*x^2/2 + Fp*x

x = 2*R;
Mo3 = vpa(subs(Mo3),3)

MoA = 200e3 == Mo3

Fp = vpa(solve(MoA,Fp),3)

%Mo1 = subs(Mo1);
%Mo2 = subs(Mo2);
%Mo3 = subs(Mo3);

%fplot(Mo1,[0,pi/2],'color','red')
%hold on
%fplot(Mo2,[0,pi/2],'color','blue')
%fplot(Mo3,[0,2*R],'color','black')
