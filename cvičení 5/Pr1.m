clear all
clc
close all
%navrh sendvicoveho nosniku zatizeneho trojuh. zat..
syms Q RA RB RC q0 m a q x Ef Iy d bet Gc S C Nc Tc Mc D f z wMo1(x) wT1(x) wMo2(x) wT2(x) X1 X2 X3 X4 X5 X6 C1 C2
%spoj zatizeni
%Q=m*10
 q0=2*Q/a;
%q0=8
 q=q0*x/a;
%q=q0
%Q=8*1800;
%rovnice SR
Fy=RA+RB+RC==Q;
MA=RB*a+RC*2*a==Q*5*a/3;
%MA=RB*a+RC*2*a==Q*3*a/2
% je nutne vyjadrit reakce jako funkce RC, protoze ta vystupuje v def. podm
RB=solve(MA,RB);
Fy=subs(Fy);
RA=solve(Fy,RA);

%rovnice vnitrnich ucinku. 
T1=RA;
Mo1=RA*x;

 T2=RA+RB-q*(x)/2;
%T2=RA+RB-q*(x-a)
Mo2=int(T2,x)+C;
r1=subs(Mo2,x,a)==0;
C=solve(r1,C);
Mo2=subs(Mo2)

%derivace
dMo1=diff(Mo1,RC);
dMo2=diff(Mo2,RC);

dT1=diff(T1,RC);
dT2=diff(T2,RC);

% def. podminka
wc=0==1/(D)*(int(Mo1*dMo1,x,0,a)+int(Mo2*dMo2,x,0,a))+bet/(Gc*S)*(int(T1*dT1,x,0,a)+int(T2*dT2,x,0,a));

Ef=200000;
Ec=250;
muc=0.06;
bet=1;
Gc=Ec/(2*(1+muc))
b=40;
%pozor h je polovina vysky jadra
h=41;
f=1;
D=Ef*2*b*f*h^2;
S=(2*h+2*f)*b;

Q=1500;
Ref=300;
RmC=4;
a=1800;

wc=subs(wc)
RC=vpa(solve(wc,RC),4)
RA=subs(RA)
RB=subs(RB)
%nyni dosadime zname reakce zpet do rovnic vnitrni ucinku pro momenty a
%pos. silu
Mo1=subs(Mo1);
Mo2=subs(Mo2);


T1=subs(T1);
T2=subs(T2);

%napeti v pasnici je dano jen napetim od ohybovho momentu, takze najdeme
%extrem
fplot(Mo1,[0 a],'color','red')
hold on
fplot(Mo2,[0 a],'color','blue')
hold off
%extrem je v druhém úseku
r3=diff(Mo2,x)==0
xmax=solve(r3,x)
%platna je jenom kladna hodnota
xmax=xmax(2,1)
Momax=abs(vpa(subs(Mo2,x,xmax),4))

sigmaFmax=Momax/(2*b*f*h)
%bezpecnost pasnic vuci MSP
kkF=Ref/sigmaFmax

%v jadre je kombinace normaloveho napeti a smyku
sigmaC1=Ec/Ef*Mo1*z/(2*b*f*h^2)
sigmaC2=Ec/Ef*Mo2*z/(2*b*f*h^2)

tauC1=T1/(2*D)*Ec*(h^2-z^2)+T1/(2*D)*Ef*(2*h*f+f^2)
tauC2=T2/(2*D)*Ec*(h^2-z^2)+T2/(2*D)*Ef*(2*h*f+f^2)

sigmaCred1=sqrt(sigmaC1^2+3*tauC1^2)
sigmaCred2=sqrt(sigmaC2^2+3*tauC2^2)

% redukovane napeti je funkce dvou promennych x a z, takze vykreslujeme ve
% 3d
fsurf(sigmaCred1,[0 a -h h])
%pro predhlednost zobrazujeme red napeti zvlast v kazdem intervalu
fsurf(sigmaCred2,[0 a -h h])

%z graf; je zrejme ze max. red. napeti v jadere je v jeho stredu pod
%podporou C
sigmaCmax=subs(sigmaCred2,z,0)
sigmaCmax=vpa(subs(sigmaCmax,x,a),4)

%bezpecnost vuci selhani jadra je tedy
kkC=RmC/sigmaCmax

%take je treba zkontrolvat tlačenou pasnici na ztratu stability
sigmaLoS=sqrt(2*f*Ef*Ec/(6*h*(1-muc^2)))
kkLoS=sigmaLoS/sigmaFmax
