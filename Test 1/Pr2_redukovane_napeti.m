clc;
clear all;
close all;

syms F Rc x a z z1 y

% urcitost i = 3-2-1 = 0 S.N

% S.R.
Rbx = F;
rc1 = 0 == Rc*3*a-F*2*a+F*a;
Rc = solve(rc1,Rc)

% usek 1 x<0,a>
N1 = F;
T1 = -Rc;
Mo1 = Rc*x

% usek 2 x<0,2*a/2>
N2 = F;
T2 = -Rc + F;
Mo2 = Rc*(a+x) - F*x;

% usek 3 x<0,a>
N3 = F;
T3 = -Rc+F-F;
Mo3 = Rc*(2*a+x) - F*(a+x) + F*x;

a = 200;
b = 80;
h = 160;
t = 20;
Iyc = 1.11e7;
zT = 60.9;
F = 50e3;
S = t*b+t*(h-t)

N1 = subs(N1);
N2 = subs(N2);
N3 = subs(N3);

T1 = subs(T1);
T2 = subs(T2);
T3 = subs(T3);

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);
Mo3 = subs(Mo3);

% vykresleni momenty pro kontrolu
%fplot(Mo1,[0,a],'color','red')
%hold on
%fplot(Mo2,[0,a],'color','blue')
%fplot(Mo3,[0,a],'color','green')


MomaxA = subs(Mo2,x,a/2)
% ohyb v bode A je 0 pusobi tam jen tah a posouvajici sila 


NA = N2
sigmaNa = vpa(NA/S) 



U1 = int(int(z1,[z,zT]),y,[-b/2,b/2])
U1max = subs(U1,z,zT-t)
U2 = int(int(z1,[z,zT-t]),y,[-t/2,t/2])
U2max = subs(U2,z1,zT-h)

%Tau = T*U/(Iyc*b)
Tau1 = T2*U1/(Iyc*b)
Tau2 = T2*(U1max+U2)/(Iyc*t)

fplot(Tau1,[zT-t,zT],'color','red')
hold on
fplot(Tau2,[zT-h,zT-t],'color','blue')

% max tau od posouvajici sily je v tezisti
Taumax = vpa(subs(Tau2,z,0),3)

% redukovane napeti podle guesta
sigmaRed = vpa(sqrt(sigmaNa^2 + 4*(Taumax)^2),3)

