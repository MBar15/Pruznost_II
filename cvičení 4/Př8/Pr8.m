clc;
clear all;
close all;

%{
    4 useky
    i = 3-2-1 = 0 S.N
    vnitrne 3 S.N.
%}

syms Tc Nc Mc E Iy x1 x phi C1 C2 C3 C4 q R phi1 d kk

Ra = 3/2*q*R;
Rb = 0;
Rbx = -1/2*q*R;

N1 = Nc*cos(phi1) + Tc*sin(phi1);
T1 = Tc*cos(phi) - Nc*sin(phi);
Mo1 = int(T1*R,phi,0,phi1) + C1;
C1 = Mc;
Mo1 = subs(Mo1);

N2 = Tc*cos(phi1) - Nc*sin(phi1) - Ra*sin(phi1)
T2 = -Tc*sin(phi) - Nc*cos(phi) - Ra*cos(phi)
Mo2 = int(T2*R,phi,0,phi1) + C2;
rdce2 = subs(Mo1,phi1,pi/2) == subs(Mo2,phi1,0);
C2 = solve(rdce2,C2)
Mo2 = subs(Mo2)

N3 = -Tc;
T3 = Rbx + Ra + Nc;
Mo3 = int(T3,x,0,x1) + C3;
rdce3 = subs(Mo3,x1,0) == subs(Mo2,phi1,pi/2);
C3 = solve(rdce3,C3);
Mo3 = subs(Mo3);

N4 = -Tc
T4 = Rbx + Ra + Nc - q*x
Mo4 = int(T4,x,0,x1) + C4;
rdce4 = subs(Mo4,x1,0) == subs(Mo3,x1,R);
C4 = solve(rdce4,C4);
Mo4 = subs(Mo4);

dMo1Nc = diff(Mo1,Nc);
dMo2Nc = diff(Mo2,Nc);
dMo3Nc = diff(Mo3,Nc);
dMo4Nc = diff(Mo4,Nc);

dMo1Tc = diff(Mo1,Tc);
dMo2Tc = diff(Mo2,Tc);
dMo3Tc = diff(Mo3,Tc);
dMo4Tc = diff(Mo4,Tc);

dMo1Mc = diff(Mo1,Mc);
dMo2Mc = diff(Mo2,Mc);
dMo3Mc = diff(Mo3,Mc);
dMo4Mc = diff(Mo4,Mc);


wNc = 0  == 1/(E*Iy) * (int(Mo1*dMo1Nc*R,phi1,0,pi/2) + int(Mo2*dMo2Nc*R,phi1,0,pi/2) + int(Mo3*dMo3Nc,x1,0,R) + int(Mo4*dMo4Nc,x1,0,R));
wTc = 0  == 1/(E*Iy) * (int(Mo1*dMo1Tc*R,phi1,0,pi/2) + int(Mo2*dMo2Tc*R,phi1,0,pi/2) + int(Mo3*dMo3Tc,x1,0,R) + int(Mo4*dMo4Tc,x1,0,R));
fiMc = 0 == 1/(E*Iy) * (int(Mo1*dMo1Mc*R,phi1,0,pi/2) + int(Mo2*dMo2Mc*R,phi1,0,pi/2) + int(Mo3*dMo3Mc,x1,0,R) + int(Mo4*dMo4Mc,x1,0,R));

Re = 300;
E = 2.1e5;
q = 15;
R = 800;
Iy=pi*d^4/64;
Wo=Iy/(d/2);
S=pi*d^2/4;
kk = 1.8;

[Q Res] = equationsToMatrix([wNc,wTc,fiMc], [Nc,Tc,Mc]);

Q = subs(Q);
Res = subs(Res);
sol = linsolve(Q,Res);
Nc = vpa(sol(1,1),3)
Tc = vpa(sol(2,1),3)
Mc = vpa(sol(3,1),3)

Mo1=subs(Mo1);
Mo2=subs(Mo2);
Mo3=subs(Mo3);
Mo4=subs(Mo4);

N1=subs(N1);
N2=subs(N2);
N3=subs(N3);
N4=subs(N4);


%pro kontrolu (odhaleni zjevnych chyb) vykreslime prubehy momentu
fplot(Mo1,[0 pi/2],'color','red')
hold on
fplot(Mo2,[0 pi/2],'color','blue')

fplot(Mo3,[0 800],'color','green')

fplot(Mo4,[0 800],'color','black')
hold off

% a normalovych sil
fplot(N1,[0 pi/2],'color','red')
hold on
fplot(N2,[0 pi/2],'color','blue')

fplot(N3,[0 800],'color','green')

fplot(N4,[0 800],'color','black')
hold off


sigmared1=abs(N1)/S+abs(Mo1)/Wo
sigmared2=abs(N2)/S+abs(Mo2)/Wo
sigmared3=abs(N3)/S+abs(Mo3)/Wo
sigmared4=abs(N4)/S+abs(Mo4)/Wo

% vykreslime je pro libovolnou velikost prumeru. jde nam o ziskani mista extremu
% ne o skutecnou velikost tohoto extremu
fplot(subs(sigmared1,d,20),[0 pi/2],'color','red')
hold on
fplot(subs(sigmared2,d,20),[0 pi/2],'color','blue')
% pro nazornost je lepsi si script zastavit zde a podivat se zvlast na
% zakrivenou cast a rovnou cast
fplot(subs(sigmared3,d,20),[0 800],'color','green')

fplot(subs(sigmared4,d,20),[0 800],'color','black')

d = vpa(solve(Re/kk==subs(sigmared1,phi1,pi/2),d),3)

