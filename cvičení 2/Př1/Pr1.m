
clc;
clear all;
close all;

syms kk Re delta b h l q alfa

Moymax = q * cosd(delta) * l^2 /2;
Mozmax = q * sind(delta) * l^2 /2;

Iy = 1/12 * b^3 * h;
Iz = 1/12 * b * h^3;

sigmay = Mozmax/Iy * b/2;
sigmaz = Moymax/Iz * h/2;

sigma = alfa*(sigmay+sigmaz) == Re/kk;

kk = 2;
delta = 45;
b = 20;
h = 30;
l = 500;
alfa = 2;
Re = 2/3 * 600;

sigma = subs(sigma);

q = abs(vpa(solve(sigma,q),3))
