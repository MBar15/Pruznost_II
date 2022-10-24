clear all;
clc;
close;

syms F R d Re kk x Rax phi Ra E Iy

F1 = F;
F2 = F/3;
a = R*sin(phi);
b = R*(1-cos(phi));

Mo1 = F1*x;
dMo1Rax = diff(Mo1,Rax);
dMo1Ra = diff(Mo1,Ra);

Mo2 = F1*(R+a) - Rax*a + Ra*(R-b);
dMo2Rax = diff(Mo2,Rax);
dMo2Ra = diff(Mo2,Ra);

Mo3 = F1*(R+b) - Rax*b + Ra*(R+a) - F2*(R-b);
dMo3Rax = diff(Mo3,Rax);
dMo3Ra = diff(Mo3,Ra);

Mo4 = F1*(R-x) + Rax*x + Ra*2*R - F2*(R-x);
dMo4Rax = diff(Mo4,Rax);
dMo4Ra = diff(Mo4,Ra);

uA = 0 == 1/(E * Iy) * (int(Mo1*dMo1Rax,x,0,R) + int(Mo2*dMo2Rax,phi,0,pi/2) + int(Mo3*dMo3Rax,phi,0,pi/2) + int(Mo4*dMo4Rax,x,0,R));
wA = 0 == 1/(E * Iy) * (int(Mo1*dMo1Ra,x,0,R) + int(Mo2*dMo2Ra*R,phi,0,pi/2) + int(Mo3*dMo3Ra*R,phi,0,pi/2) + int(Mo4*dMo4Ra,x,0,R));

Re = 250;
kk = 2;
E = 2.1e5;
d = 30;
Iy = pi*d^4/64;
R = 480;

uA = subs(uA);
wA = subs(wA);


[Q Res] = equationsToMatrix([uA,wA],[Rax,Ra]);
sols = linsolve(Q,Res);
Rax = vpa(sols(1,1),3)
Ra = vpa(sols(2,1),3)


Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);
Mo4 = subs(Mo4);

fplot(subs(Mo1,F,1),[0 R],'color', 'red')
hold on
fplot(subs(Mo2,F,1),[0 pi/2],'color', 'blue')
hold on
fplot(subs(Mo3,F,1),[0 pi/2],'color', 'green')

fplot(subs(Mo4,F,1),[0 R],'color', 'black')

%Momax = subs(Mo2, phi,0) + subs(Mo1,x,R);
%sigma = Momax*32/(pi*d^2) == Re / kk
%F=vpa(solve(sigma,F),5)


