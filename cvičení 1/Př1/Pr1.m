clc;
clear all;
close all;

% definovat symbolicke proměnné
syms F c v a sigmak Fd x E Iy Momax Wo Md





% ----------------- řesení velikosti průhybu ----------------- wFd = 3.72 > 0.5
%{

% Momenty úseků
Mo1 = -Fd*x;
dMo1 = diff(Mo1,Fd);

Mo2 = -Fd*(c+x) - F*x;
dMo2 = diff(Mo2,Fd);

% rovnice průhybu
wFd = 1/(E*Iy)*(int(Mo1*dMo1,x,0,c)+int(Mo2*dMo2,x,0,c))

% dosazení do symbolických proměných
F = 100;
c = 500;
a = 20;
v = 0.5;
Iy = 1/12*a^4;
E = 2.1e5;
Fd = 0;

% substituce hodnot místo symbolických proměnných
wFd = vpa(subs(wFd),3) 
%}

% ----------------- řesení velikosti Fd -----------------Fd = 27,1 ; Momax = 23 N*m
%{

% Momenty úseků
Mo1 = Fd*x;
dMo1 = diff(Mo1,Fd);

Mo2 = Fd*(c+x) - F*x;
dMo2 = diff(Mo2,Fd);

% rovnice průhybu
wFd = -0.5 == 1/(E*Iy)*(int(Mo1*dMo1,x,0,c)+int(Mo2*dMo2,x,0,c));

% dosazení do symbolických proměných
F = 100;
c = 500;
a = 20;
v = 0.5;
Iy = 1/12*a^4;
E = 2.1e5;

% substituce hodnot místo symbolických proměnných a řešení rovnice (solve(f))
wFd = subs(wFd);
Fd = vpa(solve(wFd,Fd),3)

%  výkreslení průběhů (zatím neznám lepší způsob)
Mo1 = subs(Mo1);
Mo2 = subs(Mo2);

fplot(Mo1,[0, c], 'color', 'red')
hold on
fplot(Mo2,[0, c], 'color', 'blue')

%}


% ----------------- výpočet napětí ----------------- % Sigma 17,2 k = 6,95
%{
sigmakk = 120;
Momax = 23e3;
a = 20;
Wo = 1/6*a^3;

sigma = Momax/ Wo

k = sigmakk/sigma
%}

% ----------------- výpočet natočení ----------------- % phiMd = 0,215 stupňů
%{

% Momenty úseků
Mo1 = Fd*x + Md;
dMo1 = diff(Mo1,Md);

Mo2 = Fd*(c+x) - F*x + Md;
dMo2 = diff(Mo2,Md);

% rovnice natočení v místě Md
phiMd = 1/(E*Iy) * (int(Mo1*dMo1,x,0,c) + int(Mo2*dMo2,x,0,c));

% dosazení do symbolických proměných
F = 100;
c = 500;
a = 20;
Iy = 1/12*a^4;
E = 2.1e5;
Fd = 27.1;
Md = 0;

% substituce hodnot místo symbolických proměnných a převod na stupně
phiMd = vpa(rad2deg(subs(phiMd)),3)
%}
