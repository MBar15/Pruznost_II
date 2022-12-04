clc;
clear all;
close all;

%{
    i = 3 - (3) = 0 S.N.
    1 průběh
%}

syms q Mk a x Iy alfa

% x <0,a>
Mo1 = q*cosd(alfa)*x^2/2;
N = q*sind(alfa)*x
Mk1 = -Mk

Re = 2/3 * 500;
q = 7;
a = 600;
kk = 1.8;
d = 60;
alfa = 30;

Mo1 = subs(Mo1)
N = subs(N)

fplot(Mo1,[0,a],'color','red')

% max zatížení bude u vetknutí
Tau = Mk1*16/(pi*d^3)
sigmaNorm = abs(subs(Mo1,x,a))*32/(pi*d^3) + subs(N,x,a)/(pi*d^2/4)

sigmaN= subs(N,x,a)/(pi*d^2/4)
sigmaM = abs(subs(Mo1,x,a))*32/(pi*d^3)

sigma = vpa(sqrt(sigmaNorm^2+4*Tau^2),3)

Mk = vpa(solve(kk == Re/sigma,Mk),3)

sigma = subs(sigma)

% Mk = 3.72e6



