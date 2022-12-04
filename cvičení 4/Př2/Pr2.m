clc;
clear all;
close all;

syms F a M1 Re E d x deltaT alfa S

Mo1 = F*x;
Mo2 = F*(a+x);

Mk1 = 0;
Mk2 = M1;
N = deltaT*alfa*E*S;

a = 300;
F = 500;
M1 = 500e3;
Re = 300;
E = 2.1e5;
d = 40;
Wo = pi*(d^3)/32;
Wk  = pi*(d^3)/16;
alfa = 12e-6; % v kelvinech na metr
S = pi*d^2/4;
kk = 1.6;


Mo1 = subs(Mo1);
Mo2 = subs(Mo2)
Mk1 = subs(Mk1);
Mk2 = subs(Mk2);



sigmaNormMax = subs(Mo2,x,a)*32/(pi*d^3) + E*alfa*deltaT;
Taumax = Mk2*16/(pi*d^3);

sigmaRed = sqrt(sigmaNormMax^2 + 4*Taumax^2)

deltaTsolved = vpa(solve(sigmaRed == Re/kk),3)

sigmaRed1 = vpa(subs(sigmaRed,deltaT,48.4))
sigmaRed2 = vpa(subs(sigmaRed,deltaT,-86.3))

% změna teploty <-48,4; 48,4>

% nutno podotknout znaménka u momentu jsou naopak oproti konvenčnímu značení,
% a sigma normmax tah-tlak- je záporné ale ono to vyjde že se to sečtě tlak
% tlak tak bych musel dát abs() je to nakonec dobře ( oba mají stejné
% znaménko) ale měly by to být trochu jinak zapsané aby to bylo správně 'gramaticky'

