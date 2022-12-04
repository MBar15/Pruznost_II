clear all
clc
close all

% vypsat vsechny konstanty, ktere budou v reseni vystupovat. Prov�d�me
% itera�n�. Nesm�me zahrnout integra�n� konstanty C1..C4
syms Q1(r) Q11(r) Q12(r) r p T D Q21(r) Q22(r) fi E mu R2 h w(r) F R1
%z�pis dif. rovnice ohybov� ��ry p�i formulaci zleva
ode = diff(Q1,r,1) ==-T/D
%rovnice popisuj�c� pr�b�hy moment�
T1=0
T2=F*R1/(2*pi*R1*r)
D=E*h^3/(12*(1-mu^2))
%vytvo�en� dvou dif. rovnic pro dva �seky
ode1=subs(ode,T,T1)
ode2=subs(ode,T,T2)

%vy�e�en� dif. rovnic. objev� se integra�n� konstanty
iQ11(r)=dsolve(ode1)
iQ12(r)=dsolve(ode2)
syms  C1 C2

ode3=1/r*diff(Q21,r,1)==iQ11(r)
ode4=1/r*diff(Q22,r,1)==iQ12(r)
iQ21(r)=dsolve(ode3)
iQ22(r)=dsolve(ode4)
syms C3 C4

rc1=iQ21(r)==r*fi
fi1=solve(rc1,fi)
Dfi1=diff(fi1,r)

rc2=iQ22(r)==r*fi
fi2=solve(rc2,fi)
Dfi2=diff(fi2,r)
% OP C3=0
C3=0

fi1=subs(fi1)
Dfi1=subs(Dfi1)
%rce napeti v usecich
sigmaR1=E*h/(2*(1-mu^2))*(Dfi1+mu*fi1/r)
sigmaT1=E*h/(2*(1-mu^2))*(fi1/r+mu*Dfi1)

sigmaR2=E*h/(2*(1-mu^2))*(Dfi2+mu*fi2/r)
sigmaT2=E*h/(2*(1-mu^2))*(fi2/r+mu*Dfi2)

rc3=subs(sigmaR1,r,R1)==subs(sigmaR2,r,R1)
rc4=subs(fi1,r,R1)==subs(fi2,r,R1)

%OP sigma R=0
rc5=0==subs(sigmaR2,r,R2)

% znamen hodnoty
F=18000
R1=100
R2=200
E=2.1e5
h=15
mu=0.3
Re=400

D=subs(D)

rc3=subs(rc3)
rc4=subs(rc4)
rc5=subs(rc5)

[Q Res]=equationsToMatrix([rc3 rc4 rc5], [C1 C2 C4])
sol=linsolve(Q,Res)
%jednotliv�m prom�nn�m p�i�ad�me p��slu�n� ��dek vektoru v�sledku. 
C1=vpa(sol(1,1),5)
C2=vpa(sol(2,1),5)
C4=vpa(sol(3,1),5)

sigmaR1=subs(sigmaR1)
sigmaT1=subs(sigmaT1)

sigmaR2=subs(sigmaR2)
sigmaT2=subs(sigmaT2)

fplot(sigmaR1,[0 R1],'color','red')
hold on

fplot(sigmaT1,[0 R1],'color','blue')

fplot(sigmaR2,[R1 R2],'color','red')


fplot(sigmaT2,[R1 R2],'color','blue')

sigmaZ=0

fplot(sigmaZ,[0 R2],'color','green')

hold off
sigmaRED=vpa(subs(sigmaR1,r,R1),5)-sigmaZ
kk=Re/sigmaRED

ode5= diff(w,r,1) ==fi1
ode6= diff(w,r,1) ==fi2
w1=dsolve(ode5)
w2=dsolve(ode6)
syms C5 C6
w1=subs(w1)
w2=subs(w2)
rc6=0==subs(w2,r,R2)
rc7=subs(w1,r,R1)==subs(w2,r,R1)

[Q Res]=equationsToMatrix([rc6 rc7], [C5 C6])
sol=linsolve(Q,Res)
%jednotliv�m prom�nn�m p�i�ad�me p��slu�n� ��dek vektoru v�sledku. 
C5=vpa(sol(1,1),5)
C6=vpa(sol(2,1),5)
w1=subs(w1)
w2=subs(w2)

fplot(w1,[0 R1],'color','green')
hold on
fplot(w2,[R1 R2],'color','green')
wmax=subs(w1,r,0)