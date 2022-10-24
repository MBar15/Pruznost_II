clear all;
clc;
close;


syms d R Re phi kk F Ra Rb Rax Rbx E Iy

Rbx = Rax;
Ra = F/2;
Rb = Ra

Mo2 = -Rax*R*sin(phi)+Ra*R*(1-cos(phi));
Mo1 = -Rax*R*sin(phi) + Rb*R*(1-cos(phi));

dMo1 = diff(Mo1,Rax);
dMo2 = diff(Mo2,Rax);

uA = 0 == 1/(E*Iy)*(int(Mo1*dMo1*R,phi,0,pi/2) + int(Mo2*dMo2*R,phi,0,pi/2));
E = 2.1e5;
d = 30;
Iy = pi*d^4/64;
R = 600;
Re = 250;
kk = 2;

uA = subs(uA);
Rax = vpa(solve(uA,Rax),3)

Mo1 = subs(Mo1)
Mo2 = subs(Mo2)

fplot(subs(Mo1,F,1),[0,pi/2],'color','red')
hold on
fplot(subs(Mo1,F,1),[0,pi/2],'color','green')

Momax = subs(Mo1,phi,pi/2)
sigma = Momax*32 / (pi*d^3) == Re/kk
F = vpa(solve(sigma,F),4)

