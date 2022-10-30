clc;
clear all;
close all;

syms P r t Re kk

%sigmared1 = P*r/t - (-P) == Re/kk

Re = 2/3*520;
t = 5;
kk = 1.5;
R = 150;
l = 12000;
r = 150;

%sigmared1 = subs(sigmared);

%P = vpa(solve(sigmared,P),3) == 7.46

sigmared2 = P*r/(2*t) - (-P) == Re/kk 

P = vpa(solve(sigmared2,P),3)





