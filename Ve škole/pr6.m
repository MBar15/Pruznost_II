clear all
clc
close all

syms F Rb Rbx R fi Iy E F2 x
d=30;
R=480;
E=210e3;
Re=250;
kk=2;
Iy=pi*d^4/64;

a=R*sin(fi);
b=R*(1-cos(fi));
Mo1=F*x;
Mo2=F*(R+a)-Rbx*a+Rb*b;
Mo3=F*(2*R-b)-F2*b+Rb*(R+a)-Rbx*(R-b);
Mo4=F*(R-x)+Rb*2*R-F2*(R+x)-Rbx*x;

dMo1Rbx=diff(Mo1,Rbx);
dMo2Rbx=diff(Mo2,Rbx);
dMo3Rbx=diff(Mo3,Rbx);
dMo4Rbx=diff(Mo4,Rbx);

dMo1Rb=diff(Mo1,Rb);
dMo2Rb=diff(Mo2,Rb);
dMo3Rb=diff(Mo3,Rb);
dMo4Rb=diff(Mo4,Rb);

ub=0==1/(E*Iy)*(int(Mo1*dMo1Rbx,x,0,R)+int(Mo2*dMo2Rbx*R,fi,0,pi/2)+int(Mo3*dMo3Rbx*R,fi,0,pi/2)+int(Mo4*dMo4Rbx,x,0,R));
vb=0==1/(E*Iy)*(int(Mo1*dMo1Rb,x,0,R)+int(Mo2*dMo2Rb*R,fi,0,pi/2)+int(Mo3*dMo3Rb*R,fi,0,pi/2)+int(Mo4*dMo4Rb,x,0,R));
F2=F/3;
ub=subs(ub)
vb=subs(vb)
[Q Res]=equationsToMatrix([ub vb], [Rbx,Rb]);

sol=linsolve(Q,Res);
Rbx=vpa(sol(1,1),5)
Rb=vpa(sol(2,1),5)

Mo1=subs(Mo1);
Mo2=subs(Mo2);
Mo3=subs(Mo3);
Mo4=subs(Mo4);
fplot(subs(Mo1,F,1),[0 R],'color','red');
hold on
fplot(subs(Mo2,F,1),[0 pi/2],'color','blue');
fplot(subs(Mo3,F,1),[0 pi/2],'color','green');
fplot(subs(Mo4,F,1),[0 R],'color','black');
hold off
Momax=abs(subs(Mo4,x,R));
sigma=Momax*32/(pi*d^3)==Re/kk;
F=vpa(solve(sigma,F),5);
Mo1=subs(Mo1);
Mo2=subs(Mo2);
Mo3=subs(Mo3);
Mo4=subs(Mo4);

%fplot(Mo1,[0 R],'color','red');
%hold on
%fplot(Mo2,[0 pi/2],'color','blue');
%fplot(Mo3,[0 pi/2],'color','green');
%fplot(Mo4,[0 R],'color','black');