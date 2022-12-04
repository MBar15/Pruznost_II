clc;
clear all;
close all;

%{
    Je to pr. 5
    i = 3 - (2+2)= 1 S.N.
    1 úsek
    Posouvající síla a její napětí Tau je zanedbána
%}

syms q a alfa Re h b E Iy Ra x Rax S Rb Rbx

Ra = q*cosd(alfa)*800/2,3
Rax = q*sind(alfa)*800/2

Mo = Ra*x - q*cosd(alfa)*x^2/2;
N = q*sind(alfa)*x - Rax;

a = 800;
b = 30;
h = 40;
alfa = 30;
Re = 250;
E = 2.1e5;
Iy = 1/12*b*h^3
q = 10;
S = b*h;

Mo = subs(Mo)
N = subs(N)

% max ohyb uprostrřed  + tahová složka
fplot(Mo,[0,a],'color','red');
hold on
fplot(N,[0,a],'color','blue');

sigmaNormMax = Mo*6/(b*h^2) + N/S,3
xmax = solve(diff(sigmaNormMax,x) == 0)
sigmaNormMax = subs(sigmaNormMax,x,xmax);

% výpočet max napětí + bezpečnost
sigmaNorm = vpa(abs(subs(Mo,x,a/2))*6/(b*h^2) + subs(N,x,a/2)/S,3)
sigmaO = vpa(abs(subs(Mo,x,a/2)*6)/(b*h^2),3)
sigmaN = subs(N,x,a/2)/S

kk = vpa(Re/sigmaNormMax,3)
% bezpečnosts je 0,28, zadáno 100*800 = 80 000 N na průřez 30*40, takže asi
% blbě zadáno pro q = 10 je bezpečnost 2.88 což dává větší smysl

