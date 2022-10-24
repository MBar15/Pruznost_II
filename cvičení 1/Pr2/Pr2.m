clc;
clear all;
close all;
 
syms kk Re E Iy d R x Ra Rb Md a q

Rb = (q*8*a^2-Md-Ra*4*a)/(3*a);

Mo1 = Ra*x - q*x^2/2 + Md;
dMo1Md = diff(Mo1,Md);

Mo2 = Ra*(a+x) + Rb*x - q*(a+x)^2/2 + Md;
dMo2Md = diff(Mo2,Md);

fiMd = 10e-3 == 1/(E*Iy) * (int(Mo1*dMo1Md,x,0,a) + int(Mo2*dMo2Md,x,0,3*a))

E = 2.1e5;
d = 30;
a = 400;
Iy = pi*d^4/64;
Md = 0;

fiMd = subs(fiMd)
Ra = solve(fiMd,Ra)

Mo1 = subs(Mo1)
Mo2 = subs(Mo2)



fplot(subs(Mo1,q,1),[0,a],'color','red')
hold on
fplot(subs(Mo2,q,1),[0, a],'color','blue')

Momax = subs(Mo2,q,1);
dMomax= 0 == diff(Momax,x)
maxX = vpa(solve(dMomax,x),3)

%  max Mo u vrubu je 166 849 N*mm
% max Mo je na x= 591
