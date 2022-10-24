clc;
clear all;
close all;

%{
    i = 3 - 3 - 1 = 1 S.N.
    2 useky
   nevim co znamena to znaceni sily F2 predpokladam ze je to ode mne
%}

syms F1 F2 Ra phi R Iy E alfa

a = R*sin(phi);
b = R-R*cos(phi)

Mo1y = F1*sind(alfa)*a + F1*cosd(alfa)*b 
Mo1z = F2*b 
Mk1 = F2*a
Mo2y = F1*sind(alfa)*(R-b) + F1*cosd(alfa)*(R+a) - Ra*a 
Mo2z = F2*(R+a) 
Mk2 = F2*(R-b)

dMo1yRa = diff(Mo1y,Ra)
dMo2yRa = diff(Mo2y,Ra)
wA = 0 == 1/(E*Iy)*(int(Mo1y*dMo1yRa,phi,0,pi/2) + int(Mo2y*dMo2yRa,phi,0,pi/2))

E = 70e3;
Re = 280;
F1 = 2500;
F2 = 800;
alfa = 15;
R = 300;
d = 35;
Iy = pi*d^4/64;

wA = subs(wA);
Ra = vpa(solve(wA,Ra))

Mo1y = subs(Mo1y);
Mo1z = subs(Mo1z);
Mo2y = subs(Mo2y);
Mo2z = subs(Mo2z);
Mk1 = subs(Mk1);
Mk2 = subs(Mk2);

sigmaNorm1 = Mo1y*32/(pi*d^3) + Mo1z*32/(pi*d^3)
sigmaNorm2 = Mo2y*32/(pi*d^3) + Mo2z*32/(pi*d^3)
Tau1 = Mk1*16/(pi*d^3)
Tau2 = Mk2*16/(pi*d^3)

sigmaRed1 = sqrt(sigmaNorm1^2+4*Tau1^2)
sigmaRed2 = sqrt(sigmaNorm2^2+4*Tau2^2)

fplot(sigmaRed1,[0,pi/2],'color','red')
hold on
fplot(sigmaRed2,[0,pi/2],'color','blue')

kk = vpa(Re/subs(sigmaRed2,phi,0),3)

% kombinované zase nevím jak sčítat

