clc;
clear all;
close all;

syms kk Re E Iy x Ra a q Ma

%Ma = q*a^2/2 - Ra * 2*a;
Ra = (q*a^2/2 - Ma)/(2*a);

Mo1 = Ma + Ra*x;
dMo1Ma = diff(Mo1,Ma);


Mo2 = Ma + Ra*(a+x) - q*x^2/2;
dMo2Ma = diff(Mo2,Ma);



fiMa = 1e-3 == 1/(E*Iy) * (int(Mo1*dMo1Ma,x,0,a) + int(Mo2*dMo2Ma,x,0,a));


kk = 2;
a = 800;
b = 30;
h = 60;
Iy = 1/12*b*h^3;
E = 2.1e5;
Re = 2/3*520;


fiMa = subs(fiMa);

Ma = solve(fiMa,Ma)
Ra = subs(Ra)

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);

%fplot(subs(Mo1,q,1),[0,a],'color','red')
%hold on 
%fplot(subs(Mo2,q,1),[0,a],'color','blue')

xmax = 0 == diff(subs(Mo2,q,1),x)
xmax = vpa(solve(xmax,x),3)
Momax = abs(subs(Mo2,x,xmax));

sigma = Momax/(1/6*b*h^2) == Re/kk;

q = vpa(solve(sigma,q),3)

% q = 24,8 
Ma = vpa(subs(Ma),3)
Ra = vpa(subs(Ra),3)
