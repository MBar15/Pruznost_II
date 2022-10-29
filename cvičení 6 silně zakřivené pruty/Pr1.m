clc;
clear all;
close all;

syms F phi R ro z y




ri = 60;
F = 50e3;
t = 20;
h = 160;
b = 80;
Re = 300

rT = (b*t*t/2 + t*(h-t)*((h-t)/2 + t)) / (b*t + t*(h-t)) +ri
S = b*t + t*(h-t)
r = vpa(S/(int(int(1/ro,ro,ri,ri+t),y,-b/2,b/2) + int(int(1/ro,ro,ri+t,ri+h),y,-t/2,t/2)),3)

e = rT - r

N = F*sin(phi);
T = -F*cos(phi);
Mo = F*rT*sin(phi);

Nmax = subs(N,phi,pi/2);
Momax = vpa(subs(Mo,phi,pi/2),3)

sigmaN = Nmax/S;
sigmaO = Momax*z/(S*e*(r-z))

z1 = r-ri-h
z2 = r-ri
fplot(sigmaO,[-118.2578, 41.74],'color','red')
hold on

sigmared = sigmaN + sigmaO;
fplot(sigmared,[-118.2578, 41.74],'color','blue')
sigmaredMax = subs(sigmared,z,41.7)
sigmaredMax1 = subs(sigmared,z,99)

kk = Re/sigmaredMax
