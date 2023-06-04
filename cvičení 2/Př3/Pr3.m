clc;
clear all;
close all;

syms F Ra Rb Rax Rbx d kk Re delta R E Iy Ra Rb Rax Rbx  

Ra = F/2;
Rb = Ra;
Rax = Rbx;

Mo1 = -Rbx*sin(delta)*R + Rb*(R-R*cos(delta));
dMo1 = diff(Mo1,Rbx);

Mo2 = -Rax*R*sin(delta) + Ra*(R - R*cos(delta));
dMo2 = diff(Mo2,Rbx);

wB = 0 == 1/(E*Iy)*(int(Mo1*dMo1*R,delta,0,pi/2) + int(Mo2*dMo2*R,delta,0,pi/2))

E = 2.1e5;
R = 600;
d = 30;
Iy = pi*d^4/64;
Re = 250;
kk = 2;

wB = subs(wB);
Rbx = vpa(solve(wB,Rbx),3)

Mo1=subs(Mo1);
Mo2=subs(Mo2);

fplot(subs(Mo1,F,1),[0 pi/2],'color','red')
hold on
fplot(subs(Mo2,F,1),[0 pi/2],'color','green')

Momax = subs(Mo1,delta,pi/2);

sigma = Momax*32/(pi*d^3) == Re/kk;

F = vpa(solve(sigma,F),4)
