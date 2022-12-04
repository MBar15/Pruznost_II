close all
clear all
clc

syms F S1 S2 DT DN phi G m g  r v h a c d e b




a = F*cos(phi) - G*sind(phi)
v = sqrt((2*r/m)*(F*sind(phi) + G*cosd(phi) - G))
G = m*g;
g = 9.81;
F = 555
r = 0.35;
m = 75
a = 0.45;
b = 0.25;
c = 0.1;
d = 0.2;
e = 0.15;


v= subs(v)

%fplot(v,[0, 74.0568],'color','red')
%xlabel('Úhel [°]')
%ylabel('v [m/s]')


DT = m*a;
DN = m*v^2/r

%x = 0 == F-S1*sind(phi) - S2*sind(phi) - DT*cosd(phi) + DN*sind(phi)
%y = 0 == S1*cosd(phi) + S2*cosd(phi) - G -DN*cosd(phi) + DT*sin(phi)
%Ma = 0 == F*b/2 - G*a/2 + S1*cosd(phi)*c + S2*cosd(phi)*(c+d) + DT*cosd(phi)*b/2 + DT*sind(phi)*a/2 + DN*sind(phi)*b/2 - DN*cos(phi)*a/2

x = 0 == F - S1*sind(phi) - S2*sind(phi) - DT*cosd(phi) + DN*sind(phi)
y = 0 == S1*cosd(phi) +S2*cosd(phi) - G -DN*cosd(phi) + DT*sind(phi)
Ma = 0 == F*b/2 - G*a/2 + S1*cosd(phi)*c + S2*cosd(phi)*(c+d) + DT*cosd(phi)*b/2 + DT*sind(phi)*a/2 + DN*sind(phi)*b/2 - DN*cosd(phi)*a/2


[Q Res] = equationsToMatrix([y,Ma], [S1, S2])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q,Res)
S1 = sol(1,1)
S2 = sol(2,1)

%74.0568
fplot(S1,[0,90],'color','green')
hold on
fplot(S2,[0, 90],'color','blue')
