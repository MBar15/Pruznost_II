close all
clear all
clc

syms F G beta phi r m Dn Dt ht bt s  Sh Sd 

a =  (F*sind(beta+phi) - G*cosd(phi))/m
v = sqrt((2*r)/m*(F*(cosd(beta) - cosd(beta+phi)) - G*sind(phi)))

x = 0 == -Dt + F*sind(beta+phi) - G*cosd(phi)
y = 0 == Sh + Sd - F*cosd(beta+phi) - G*sind(phi) + Dn
M = 0 == G*bt + Dt*cosd(phi)*bt + Dt*sind(phi)*ht - F*cosd(beta)*s + Sh*cosd(phi)*s + Dn*(cosd(phi)*ht - sind(phi)*bt)

m = 800;
bt = 0.4;
r = 1.8;
F = 8500;
beta = 25;
ht = 0.5;
s = 0.9;
phiH = 59.44;
g = 9.81;
G = m*g
Dn = m*v^2/r
%Dt = a/r*r*m

lul = vpa(subs(subs(Dn),phi,phiH),3)
test =vpa(subs(subs(Dt),phi,phiH),3)

[Q Res] = equationsToMatrix([x,y, M], [Sh, Sd,Dt])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q,Res)
Sh = sol(1,1)
Sd = sol(2,1)

Sh = subs(Sh)
Sd = subs(Sd)

Sh1 = vpa(subs(Sh,phi,phiH),3)
Sd1 = vpa(subs(Sd,phi,phiH),3)
Dt1 = vpa(subs(subs(sol(3,1),phi,phiH)),3)


fplot(Sh,[0,90],'color','red')
hold on
fplot(Sd,[0,90],'color','blue')
