clc;
clear all;
close all;

syms F a b h R Re Iy E Ma phi x rho z 

Mo1 = Ma + F/2*R*sin(phi);
Mo2 = Ma + F/2*R;

dMo1Ma = diff(Mo1,Ma);
dMo2Ma = diff(Mo2,Ma);

fiA = 0 == 4/(E*Iy)*(int(Mo1*dMo1Ma*R,phi,0,pi/2) + int(Mo2*dMo2Ma,x,0,a/2))

F = 50e3;
a = 100;
b = 10;
h = 100;
R = 100;
Re = 2/3*500;
E = 2.1e5;
Iy = 1/12*b*h^3;

fiA = subs(fiA);
Ma = vpa(solve(fiA,Ma),3)

Mo1 = subs(Mo1);
Mo2 = subs(Mo2);

%fplot(Mo1,[0,pi/2],'color','red')
%hold on
%fplot(Mo2,[0,a/2],'color','green')

Momax =  subs(Mo1,phi,0)

S = b*h;
Ri = R-h/2;
Rn = vpa(S/int(int(1/rho,Ri,Ri+h),-b/2,b/2),3)

e = vpa(R-Rn,3)

sigmaN = F/2/S;
sigmaMo = Momax*z/(S*e*(Rn-z))

z1 =Rn-R-h/2
z2 =Rn-Ri
fplot(sigmaMo,[-58.97,41.0239])



sigmaRedMax = subs(sigmaMo,z,41.0239)

kk = vpa(Re/abs(sigmaRedMax))

