clc;
clear all;
close all;

syms R t rhoMat rhoKap z g r P l Nm S fi sigmaM sigmaT sigmaR

%-------------- dolni pulkulove dno-----------------
% r je funkce z
rc1=1==((R-z)/R)^2+(r/R)^2;
r=solve(rc1,r);
r=r(1,1);
%meridiamova sila musi musi vyrovnat tihu kapaliny tihu mat., silu od
%hydrostat tlaku a tlak

% nezanedbam hmotnost nadoby ale zanedb8m hmotnost dna y duvodu sloziteho
% vypoctu

% tiha mat
Vmat=2*pi*r*z*t
Gmat=rhoMat*g*Vmat

%tiha kapaliny
Vkap=1/3*pi*z^2*(3*R-z)
Gkap=rhoKap*g*Vkap

% hydrostat tlak
Ph = rhoKap*g*(l+R-z)
Fh = Ph*pi*r^2

% sila od pretlaku
Fp = P*pi*r^2

% meridiamova sila
Fm = Nm*cos(fi)
Fm=subs(Fm,cos(fi),r/R)
Smat=2*pi*r*t
nap=sigmaM==Nm/Smat
nap=subs(nap)
Nm=solve(nap,Nm)
Fm=subs(Fm)

% silova rovnovaha
rc2 = Fm == Gkap + Fh + Gmat + Fp

rhoKap = 1800;
rhoMat = 7850;
R = 0.4;
l = 3;
g = 9.81;
P = 0.8e6
Re = 350e6;
kk = 1.5

rc2 = subs(rc2)

sigmaM = solve(rc2,sigmaM)
sigmaM1=subs(sigmaM,t,0.001)
fplot(sigmaM1,[0, 0.4],'color','red')

% laplaceova rovnice
rM = R;
rT = R;

pCelk = P + Ph
rc3 = pCelk/t == sigmaM/rM + sigmaT/rT
rc3 = subs(rc3)
sigmaT = solve(rc3,sigmaT)
sigmaT1=subs(sigmaT,t,0.001)
fplot(sigmaT1,[0 0.4],'color','green')

% simgaT a sigmaM je v 0 stejne a to 172 MPa
z = 0
sigmaRED1=sigmaT-(-pCelk)==Re/kk;
sigmaRED1=subs(sigmaRED1)
t1=vpa(solve(sigmaRED1,t),5)*1000 % 8 mm



% valcova cast stejně jako předešlé příklady













