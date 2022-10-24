clc;
clear all;
close all;

syms R a Re m d G fi p

Mo1 = -G*R*sind(fi);
Mo2 = -G*R;
Mo3 = -G*R*cosd(fi);

m = 120;
G = 9.81*m;
kk = 2;
R = 800;
a = 200;
Re = 400;
Wo = pi*d^3/32;

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);

%fplot(Mo1,[0,90],'color','red')
%hold on
%fplot(Mo2,'color','blue')
%fplot(Mo3,[0,90],'color','black')

Momax = abs(Mo2)

Wo = pi*d^3/32  == (Momax)/(Re/kk)
d = vpa(solve(Wo,d),3)