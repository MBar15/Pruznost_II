clc;
clear all;
close all;

syms E R Ra Rax q x Iy E fi kk sigmay Fd

% pro začátek nevím k čemu je sigmay a mý a jaký vztah je mezi nimi pro
% sigma dov, a tak je sigma dov pouze sigmay

a = R - R*cos(fi);
b = R*sin(fi);

Mo1 = Ra*a + Rax*b
dMo1 = diff(Mo1,Rax);
dMo1Fd = diff(Mo1,Fd);

Mo2 = Ra*(R + x) + Rax*R - q*x^2/2;
dMo2 = diff(Mo2,Rax);
dMo2Fd = diff(Mo2,Fd);

uA = 1/(E*Iy) * (int(Mo1*dMo1*R,fi,0,pi/2) + int(Mo2*dMo2,x,0,R));

E = 100e3;
my = 0.3;
D = 40;
Iy = pi*D^4/64;
R = 400;
q = 10;
Rax = 0;
Ra = (q*2*R*2*R/2)/(3*R);
sigmay = 770;

uA = vpa(subs(uA),3)

fplot(subs(Mo1),[0,pi/2],'color','red')
hold on
fplot(subs(Mo2),[0,R],'color','blue')



xMax = diff(Mo2,x) == 0;
xMax = vpa(solve(subs(xMax),x),3)
Momax = subs(Mo2)
Momax = subs(Momax,x,xMax)

sigma = sigmay/kk == Momax*32/(pi*D^3)
kk = vpa(solve(sigma,kk),3)