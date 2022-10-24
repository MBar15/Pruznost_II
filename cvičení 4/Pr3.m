clc;
clear all;
close all;

%{
    Je to pr. 5
    i = 3 - (2+2)= 1 S.N.
    1 úsek
    Posouvající síla a její napětí Tau je zanedbána
%}

syms q a alfa Re h b E Iy Ra x Rax S

% x <0,a>

Mo = Ra*x - q*cosd(alfa)*x^2/2;
N = q*sind(alfa)*x - Rax;

dMoRa = diff(Mo,Ra);
dNRax = diff(N,Rax)

% castiglián pro ohyb
wA = 0 == 1/(E*Iy)*int(Mo*dMoRa,x,0,a)
% castiglián pro tah-tlak
xA = 0 == 1/(E*S)*(int(N*dNRax,x,0,a))

a = 800;
b = 30;
h = 40;
alfa = 30;
Re = 250;
E = 2.1e5;
Iy = 1/12*b*h^3
q = 100;
S = b*h;

% výpočet neznámých
wA = subs(wA);
xA = subs(xA);
Ra = vpa(solve(wA,Ra),3)
Rax = vpa(solve(xA,Rax),3)

Mo = subs(Mo);
N = subs(N);

% max ohzb na konci  + tahová složka
fplot(Mo,[0,a],'color','red');
hold on
fplot(N,[0,a],'color','blue');

% výpočet max napětí + bezpečnost
sigmaNorm = vpa(abs(subs(Mo,x,a))*6/(b*h^2) + subs(N,x,a)/S,3)
sigmaO = vpa(abs(subs(Mo,x,a)*6)/(b*h^2),3)
sigmaN = vpa(subs(N,x,a)/S,3)

kk = Re/sigmaNorm

% bezpečnosts je 0,28, zadáno 100*800 = 80 000 N na průřez 30*40, takže asi
% blbě zadáno pro q = 10 je bezpečnost 2.83 což dává větší smysl

