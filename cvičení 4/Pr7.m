clc;
clear all;
close all;


%{
    i = 3-3 = 0 S.N.
    4 useky

    osa y je ve smeru F2
    osa z je ve smeru F1
    zanedbavam napeti od posouvajici slozky
    material neni zadany volim 11 600
%}


syms F1 F2 a M1 q kk x d

% x <0,a)
Mo1y = -q*x^2/2;
Mo1z = 0;
Mk1 = -M1;

% x <0,a>
Mo2y = -q*a*(a/2+x);
Mo2z = 0;
Mk2 = -M1;

% x <0,2>
Mo3y = M1+F1*x-q*a*(x);
Mo3z = 0;
Mk3 = -q*a*(a+a/2);

% x <0,a>
Mo4y = M1+F1*(x+a)-q*a*(x+a);
Mo4z = F2*x;
Mk4 = -q*a*(a+a/2);

F1 = 3500;
F2 = 5200;
a = 350;
M1 = 180e3;
q = 15;
kk = 1.8;
E = 2.1e5;
Re = 400;

Mo1y = subs(Mo1y);
Mo1z = subs(Mo1z);
Mk1 = subs(Mk1);
Mo2y = subs(Mo2y);
Mo2z = subs(Mo2z);
Mk2 = subs(Mk2);
Mo3y = subs(Mo3y);
Mo3z = subs(Mo3z);
Mk3 = subs(Mk3);
Mo4y = subs(Mo4y)
Mo4z = subs(Mo4z)
Mk4 = subs(Mk4)


% kontrola Moment; jestli to dava smysl
%Moy
%fplot(Mo1y,[0,a],'color','red')
%hold on
%fplot(Mo2y,[0,a],'color','blue')
%fplot(Mo3y,[0,a],'color','green')
%fplot(Mo4y,[0,a],'color','black')

%Moz
%fplot(Mo1z,[0,a],'color','red')
%hold on
%fplot(Mo2z,[0,a],'color','blue')
%fplot(Mo3z,[0,a],'color','green')
%fplot(Mo4z,[0,a],'color','black')

%Mk
%fplot(Mk1,[0,a],'color','red')
%hold on
%fplot(Mk2,[0,a],'color','blue')
%fplot(Mk3,[0,a],'color','green')
%fplot(Mk4,[0,a],'color','black')

sigmaNorm1 = Mo1y*32/(pi*d^3) + Mo1z*32/(pi*d^3);
sigmaNorm2 = Mo2y*32/(pi*d^3) + Mo2z*32/(pi*d^3);
sigmaNorm3 = Mo3y*32/(pi*d^3) + Mo3z*32/(pi*d^3);
sigmaNorm4 = Mo4y*32/(pi*d^3) + Mo4z*32/(pi*d^3);

Tau1 = Mk1*16/(pi*d^3);
Tau2 = Mk2*16/(pi*d^3);
Tau3 = Mk3*16/(pi*d^3);
Tau4 = Mk4*16/(pi*d^3);

sigmaRed1 = sqrt(sigmaNorm1^2 + 4*Tau1^2);
sigmaRed2 = sqrt(sigmaNorm2^2 + 4*Tau2^2);
sigmaRed3 = sqrt(sigmaNorm3^2 + 4*Tau3^2)
sigmaRed4 = sqrt(sigmaNorm4^2 + 4*Tau4^2)

fplot(subs(sigmaRed1,d,1),[0,a],'color','red')
hold on
fplot(subs(sigmaRed2,d,1),[0,a],'color','blue')
fplot(subs(sigmaRed3,d,1),[0,a],'color','green')
fplot(subs(sigmaRed4,d,1),[0,a],'color','black')

% max sigmared nevim jak zjistit secist je jen nemuzu jelikoz je to ohyb a
% krut na kulatine bud je to Mo1y+Mk + Moz(0) nebo Mo1z+Mk + Moy(0)
% nebo kombinace k*Mo1z + k*Mo1y + Mk














