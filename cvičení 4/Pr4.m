clc;
clear all;
close all;

%{
    i = 3 - (3) = 0 S.N.
    1 průběh
%}
% x <0,a>

syms q Mk a x Iy 

Mo1 = q*x^2/2;
Mk1 = -Mk

Re = 2/3 * 500;
q = 7;
a = 600;
kk = 1.8
d = 60

Mo1 = subs(Mo1);

fplot(Mo1,[0,a],'color','red')

Tau = Mk1*16/(pi*d^3)
sigmaNorm = abs(subs(Mo1,x,a))*32/(pi*d^3)

sigma = vpa(sqrt(sigmaNorm^2+4*Tau^2),3)

Mk = vpa(solve(kk == Re/sigma,Mk),3)

sigma = subs(sigma)

% Mk = 3.72e6



