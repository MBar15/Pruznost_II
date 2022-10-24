clc;
clear all;
close all;


% teziste yt = 48.1818
t = 20;
h = 160;
b = 80;

yt = ((h-t)*t*(h-t)/2 + (b*t)*t/2)/((h-t)*t+b*t)


% Kvadraticky moment
Iy = 1/12*b*t^3 + 1/12*t*(h-t)^3 + b*t*(t/2-yt)^2 + (h-t)*t*((h-t)/2+t - yt)^2


syms l Re q k Rb x E

Mo1 = Rb*x - q*x^2/2
dMo1Rb = diff(Mo1,Rb)

wB = -Rb/k == 1/(E*Iy) * int(Mo1*dMo1Rb,x,0,l);

l = 800;
Re = 350;
q = 200;
k = 10e4;
E = 2.1e5;

wB= subs(wB)
Rb = vpa(solve(wB,Rb),3)


wB = 1/(E*Iy) * int(Mo1*dMo1Rb,x,0,l); 
wB = subs(wB)

%Mo1 = subs(Mo1)
%fplot(Mo1,[0,l],'color','red');

% Momax = -22 110 800 N*mm

Momax = abs(subs(Mo1,x,l));

SigmaMaxTlak = Momax/Iy*yt;

kk = Re/SigmaMaxTlak;

