clc;
clear all;
close all;

syms kk d Re R F Rbx Rb x Iy alfa E

%b = R*sin(alfa);
%c = R- R*cos(alfa)

Mo1 = F * x;
dMo1Rb= diff(Mo1,Rb);
dMo1Rbx= diff(Mo1,Rbx);

Mo2 = F*(R+R*sin(alfa)) - Rbx*R*sin(alfa) + Rb*(R-R*cos(alfa));
dMo2Rb= diff(Mo2,Rb);
dMo2Rbx= diff(Mo2,Rbx);

Mo3 = F*(R+R*cos(alfa)) + Rb*(R+R*sin(alfa)) - Rbx*R*cos(alfa) - F/3*(R-R*cos(alfa));
dMo3Rb= diff(Mo3,Rb);
dMo3Rbx= diff(Mo3,Rbx);

Mo4 = F*(R-x) + Rb*2*R + Rbx*x - F/3*(R+x);
dMo4Rb= diff(Mo4,Rb);
dMo4Rbx= diff(Mo4,Rbx);

wA = 0 == 1/(E*Iy) * (int(Mo1*dMo1Rb,x,0,R) + int(Mo2*dMo2Rb*R,alfa,0,pi/2) + int(Mo3*dMo3Rb*R,alfa,0,pi/2) + int(Mo4*dMo4Rb,x,0,R))
uA = 0 == 1/(E*Iy) * (int(Mo1*dMo1Rbx,x,0,R) + int(Mo2*dMo2Rbx*R,alfa,0,pi/2) + int(Mo3*dMo3Rbx*R,alfa,0,pi/2) + int(Mo4*dMo4Rbx,x,0,R))

E = 2.1e5;
R = 480;
d = 30;
Iy = pi*d^4/64;
Re = 250;
kk = 2;

[Q, Res] = equationsToMatrix([wA, uA], [Rb, Rbx]);
Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q, Res);
Rb = vpa(sol(1,1));
Rbx = vpa(sol(2,1));

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);
Mo4 = subs(Mo4);


%fplot(subs(Mo1,F,1),[0,R],'color','blue');
%hold on 
%fplot(subs(Mo2,F,1),[0,pi/2],'color','green');
%fplot(subs(Mo3,F,1),[0,pi/2],'color','black');
%fplot(subs(Mo4,F,1),[0,R],'color','red');

Momax = subs(Mo1,x,R);

sigma = Momax*32/(pi*d^3) == Re/kk;
F= vpa(solve(sigma,F),3)
Mo1=subs(Mo1);
Mo2=subs(Mo2);
Mo3=subs(Mo3);
Mo4=subs(Mo4);

fplot(Mo1,[0 R],'color','red');
hold on
fplot(Mo2,[0 pi/2],'color','blue');
fplot(Mo3,[0 pi/2],'color','green');
fplot(Mo4,[0 R],'color','black');
