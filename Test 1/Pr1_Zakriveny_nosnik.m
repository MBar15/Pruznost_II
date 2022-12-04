clc;
clear all;
close all;

syms fi q Ra x a h b Iy E Re R2 R1 Rax R rho y z

%{  
    nosnik je silne zakriveny u R2,kdy R2/h = 80/120 < 1
    nosnik je slabe zakriveny u R1,kdy R1/h = 700/120 > 5
    i = 3-3-2 = 2 S.N. 
    zanedbavam vliv posouvajici sily
%}



% x<0,a>
Mo1 = -q*x^2/2 + Ra*x;

% fi<0,pi/2>
N2 = Rax*cos(fi) - Ra*sin(fi) + q*sin(fi);
T2 = -Ra*sin(fi) +q*cos(fi) + Rax*sin(fi)
Mo2 = Ra*(a+R2*sin(fi)) - q*a*(a/2+R2*sin(fi)) - Rax*(R2-R2*cos(fi))

% fi<0,pi/2>
Mo3 = Ra*(a+R2+R1*(1-cos(fi))) - q*a*(a/2+R2+R1*(1-cos(fi))) - Rax*(R2+R1*sin(fi))

% fi<0,pi/2>
Mo4 = Ra*(a+R2+R1+R1*sin(fi)) - q*a*(a/2+R2+R1+R1*sin(fi)) - Rax*(R2+R1*cos(fi))

dMo1Ra = diff(Mo1,Ra);
dMo2Ra = diff(Mo2,Ra);
dMo3Ra = diff(Mo3,Ra);
dMo4Ra = diff(Mo4,Ra);

dMo1Rax = diff(Mo1,Rax);
dMo2Rax = diff(Mo2,Rax);
dMo3Rax = diff(Mo3,Rax);
dMo4Rax = diff(Mo4,Rax);

wA = 0 == 1/(E*Iy)*(int(Mo1*dMo1Ra,x,0,a) + int(Mo2*dMo2Ra*R2,fi,0,pi/2) + int(Mo3*dMo3Ra*R1,fi,0,pi/2) + int(Mo4*dMo4Ra*R1,fi,0,pi/2))
uA = 0 == 1/(E*Iy)*(int(Mo1*dMo1Rax,x,0,a) + int(Mo2*dMo2Rax*R2,fi,0,pi/2) + int(Mo3*dMo3Rax*R1,fi,0,pi/2) + int(Mo4*dMo4Rax*R1,fi,0,pi/2))

a = 500;
b = 40;
h = 120;
R1 = 700;
R2 = 80;
t = 2;
kk = 1.5
Re = 300;
E = 2.1e5;
Iy = 1/12*b*h^3 - 1/12*(b-2*3*t)*(h-2*t)^3;
S = b*h-(b-2*3*t)*(h-2*t)

[Q Res] = equationsToMatrix([wA, uA], [Ra, Rax])
Q = subs(Q);
Res = subs(Res)
sol = linsolve(Q,Res);
Ra= sol(1,1)
Rax = sol(2,1)

Mo1 =subs(Mo1)
Mo2 =subs(Mo2)
Mo3 =subs(Mo3)
Mo4 =subs(Mo4)

%fplot(subs(Mo1,q,1),[0,a],'color','red')
%hold on
%fplot(subs(Mo2,q,1),[0,pi/2],'color','blue')
%hold on
%fplot(subs(Mo3,q,1),[0,pi/2],'color','green')
%fplot(subs(Mo4,q,1),[0,pi/2],'color','black')

Ri = R2-h/2
Ri2 = Ri+t;
Rn = vpa(S/(int(int(1/rho,rho,Ri,Ri+h),y,-b/2,b/2) - int(int(1/rho,rho,Ri2,Ri2+h-t),y,-(b-2*3*t)/2,(b-2*3*t)/2)))
e = vpa(R2 - Rn,3)

Mo2Max = subs(Mo2,fi,0)

sigma = Mo2Max*z/(S*e*(Rn-z))
sigmaN =N2/S
sigmaN = subs(sigmaN,fi,0)
sigmaN = subs(sigmaN)
z1 = Ri+h-Rn
z2 = Ri

fplot(subs(sigma,q,1),[-80,35.20])
hold on
%fplot(subs(sigmaN,q,1),[0,pi/2],'color','red')

%max je v 35.50

sigmaMax = subs(sigma,z,35.20)
rc1 = sigmaMax + sigmaN  == Re/kk
q = solve(rc1,q) % 46.89 zanedbavam vliv normalovych sil, jednat protoze jsem to nestihl a chyba by mela byt mala
