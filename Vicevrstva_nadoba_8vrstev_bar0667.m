clc;
clear all;
close all

% Osm vrstev nádoby ze stejného materiálu oce 11 600
% tlak ve vnitr nadoby p

syms r r1 r2 r3 r4 r5 r6 r7 r8 r9 A1 A2 A3 A4 A5 A6 A7 A8 B1 B2 B3 B4 B5 B6 B7 B8 p ps2 ps3 ps4 ps5 ps6 ps7 ps8 E mu Re kk delt1 delt2 delt3 delt4 delt5 delt6 delt7 

%------------ rovnice raddíálního a tečného napětí ----------------------
sigmaR1 = A1 + B1/r^2;
sigmaT1 = A1 - B1/r^2;

sigmaR2 = A2 + B2/r^2;
sigmaT2 = A2 - B2/r^2;

sigmaR3 = A3 + B3/r^2;
sigmaT3 = A3 - B3/r^2;

sigmaR4 = A4 + B4/r^2;
sigmaT4 = A4 - B4/r^2;

sigmaR5 = A5 + B5/r^2;
sigmaT5 = A5 - B5/r^2;

sigmaR6 = A6 + B6/r^2;
sigmaT6 = A6 - B6/r^2;

sigmaR7 = A7 + B7/r^2;
sigmaT7 = A7 - B7/r^2;

sigmaR8 = A8 + B8/r^2;
sigmaT8 = A8 - B8/r^2;

%------------ rovnice raddiálního napětí ------------------------
sigmaR1R1 = subs(sigmaR1,r,r1);

sigmaR1R2 = subs(sigmaR1,r,r2);
sigmaR2R2 = subs(sigmaR2,r,r2);

sigmaR2R3 = subs(sigmaR2,r,r3);
sigmaR3R3 = subs(sigmaR3,r,r3);

sigmaR3R4 = subs(sigmaR3,r,r4);
sigmaR4R4 = subs(sigmaR4,r,r4);

sigmaR4R5 = subs(sigmaR4,r,r5);
sigmaR5R5 = subs(sigmaR5,r,r5);

sigmaR5R6 = subs(sigmaR5,r,r6);
sigmaR6R6 = subs(sigmaR6,r,r6);

sigmaR6R7 = subs(sigmaR6,r,r7);
sigmaR7R7 = subs(sigmaR7,r,r7);

sigmaR7R8 = subs(sigmaR7,r,r8);
sigmaR8R8 = subs(sigmaR8,r,r8);



%------------ rovnice tečného napětí ----------------------
sigmaT1R1 = subs(sigmaT1,r,r1);

sigmaT1R2 = subs(sigmaT1,r,r2);
sigmaT2R2 = subs(sigmaT2,r,r2);

sigmaT2R3 = subs(sigmaT2,r,r3);
sigmaT3R3 = subs(sigmaT3,r,r3);

sigmaT3R4 = subs(sigmaT3,r,r4);
sigmaT4R4 = subs(sigmaT4,r,r4);

sigmaT4R5 = subs(sigmaT4,r,r5);
sigmaT5R5 = subs(sigmaT5,r,r5);

sigmaT5R6 = subs(sigmaT5,r,r6);
sigmaT6R6 = subs(sigmaT6,r,r6);

sigmaT6R7 = subs(sigmaT6,r,r7);
sigmaT7R7 = subs(sigmaT7,r,r7);

sigmaT7R8 = subs(sigmaT7,r,r8);
sigmaT8R8 = subs(sigmaT8,r,r8);



%------------ nevím jak se to nazývá ----------------------
sigmaZ1R2 = (sigmaR1R2 + sigmaT1R2)/2
sigmaZ2R2 = (sigmaR2R2 + sigmaT2R2)/2

sigmaZ2R3 = (sigmaR2R3 + sigmaT2R3)/2
sigmaZ3R3 = (sigmaR3R3 + sigmaT3R3)/2

sigmaZ3R4 = (sigmaR3R4 + sigmaT3R4)/2
sigmaZ4R4 = (sigmaR4R4 + sigmaT4R4)/2

sigmaZ4R5 = (sigmaR4R5 + sigmaT4R5)/2
sigmaZ5R5 = (sigmaR5R5 + sigmaT5R5)/2

sigmaZ5R6 = (sigmaR5R6 + sigmaT5R6)/2
sigmaZ6R6 = (sigmaR6R6 + sigmaT6R6)/2

sigmaZ6R7 = (sigmaR6R7 + sigmaT6R7)/2
sigmaZ7R7 = (sigmaR7R7 + sigmaT7R7)/2

sigmaZ7R8 = (sigmaR7R8 + sigmaT7R8)/2
sigmaZ8R8 = (sigmaR8R8 + sigmaT8R8)/2

%------------ rovnice napětí pro výpočet ----------------------

rce1 = subs(sigmaR1,r,r1) == -p
rce2 = subs(sigmaR1,r,r2) == -ps2
rce3 = subs(sigmaR2,r,r2) == -ps2
rce4 = subs(sigmaR2,r,r3) == -ps3
rce5 = subs(sigmaR3,r,r3) == -ps3
rce6 = subs(sigmaR3,r,r4) == -ps4
rce7 = subs(sigmaR4,r,r4) == -ps4
rce8 = subs(sigmaR4,r,r5) == -ps5
rce9 = subs(sigmaR5,r,r5) == -ps5
rce10 = subs(sigmaR5,r,r6) == -ps6
rce11 = subs(sigmaR6,r,r6) == -ps6
rce12 = subs(sigmaR6,r,r7) == -ps7
rce13 = subs(sigmaR7,r,r7) == -ps7
rce14 = subs(sigmaR7,r,r8) == -ps8
rce15 = subs(sigmaR8,r,r8) == -ps8
rce16 = subs(sigmaR8,r,r9) == 0

%------------ přesahy ----------------------

u12 = r2/E*(sigmaT1R2 - mu*(sigmaR1R2 + sigmaZ1R2))
u22 = r2/E*(sigmaT2R2 - mu*(sigmaR2R2 + sigmaZ2R2))

u23 = r3/E*(sigmaT2R3 - mu*(sigmaR2R3 + sigmaZ2R3))
u33 = r3/E*(sigmaT3R3 - mu*(sigmaR3R3 + sigmaZ3R3))

u34 = r4/E*(sigmaT3R4 - mu*(sigmaR3R4 + sigmaZ3R4))
u44 = r4/E*(sigmaT4R4 - mu*(sigmaR4R4 + sigmaZ4R4))

u45 = r5/E*(sigmaT4R5 - mu*(sigmaR4R5 + sigmaZ4R5))
u55 = r5/E*(sigmaT5R5 - mu*(sigmaR5R5 + sigmaZ5R5))

u56 = r6/E*(sigmaT5R6 - mu*(sigmaR5R6 + sigmaZ5R6))
u66 = r6/E*(sigmaT6R6 - mu*(sigmaR6R6 + sigmaZ6R6))

u67 = r7/E*(sigmaT6R7 - mu*(sigmaR6R7 + sigmaZ6R7))
u77 = r7/E*(sigmaT7R7 - mu*(sigmaR7R7 + sigmaZ7R7))

u78 = r8/E*(sigmaT7R8 - mu*(sigmaR7R8 + sigmaZ7R8))
u88 = r8/E*(sigmaT8R8 - mu*(sigmaR8R8 + sigmaZ8R8))


%------------ rovnice z přesahu (delty) ----------------------

rce17 = delt1 == -u12 + u22
rce18 = delt2 == -u23 + u33
rce19 = delt3 == -u34 + u44
rce20 = delt4 == -u45 + u55
rce21 = delt5 == -u56 + u66
rce22 = delt6 == -u67 + u77
rce23 = delt7 == -u78 + u88

% ------------ dovolené napětí -------------------------------

rce24 =  sigmaT1R1 - sigmaR1R1 == Re/kk
rce25 =  sigmaT2R2 - sigmaR2R2 == Re/kk
rce26 =  sigmaT3R3 - sigmaR3R3 == Re/kk
rce27 =  sigmaT4R4 - sigmaR4R4 == Re/kk
rce28 =  sigmaT5R5 - sigmaR5R5 == Re/kk
rce29 =  sigmaT6R6 - sigmaR6R6 == Re/kk
rce30 =  sigmaT7R7 - sigmaR7R7 == Re/kk
rce31 =  sigmaT8R8 - sigmaR8R8 == Re/kk


%------------- parametry -------------------------------------

r1 = 300;
r2 = 400;
r3 = 500;
r4 = 600;
r5 = 700;
r6 = 800;
r7 = 900;
r8 = 1000;
r9 = 1100;
E = 2.11e5;
Re = 3/2*600;
mu = 0.3;
kk = 1.7;

% ----------- řešení rovnic ---------------------------------

[Q Res] = equationsToMatrix([rce1 rce2 rce3 rce4 rce5 rce6 rce7 rce8 rce9 rce10 rce11 rce12 rce13 rce14 rce15 rce16 rce17 rce18 rce19 rce20 rce21 rce22 rce23 rce24 rce25 rce26 rce27 rce28 rce29 rce30 rce31], [A1 B1 A2 B2 A3 B3 A4 B4 A5 B5 A6 B6 A7 B7 A8 B8 ps2 ps3 ps4 ps5 ps6 ps7 ps8 delt1 delt2 delt3 delt4 delt5 delt6 delt7 p])
Q = subs(Q)
Res = subs(Res)
sol = linsolve(Q, Res)

% ----------- dosazení výsledků
A1 = vpa(sol(1,1),3)
B1 = vpa(sol(2,1),3)
A2 = vpa(sol(3,1),3)
B2 = vpa(sol(4,1),3)
A3 = vpa(sol(5,1),3)
B3 = vpa(sol(6,1),3)
A4 = vpa(sol(7,1),3)
B4 = vpa(sol(8,1),3)
A5 = vpa(sol(9,1),3)
B5 = vpa(sol(10,1),3)
A6 = vpa(sol(11,1),3)
B6 = vpa(sol(12,1),3)
A7 = vpa(sol(13,1),3)
B7 = vpa(sol(14,1),3)
A8 = vpa(sol(15,1),3)
B8 = vpa(sol(16,1),3)
ps2 = vpa(sol(17,1),3)
ps3 = vpa(sol(18,1),3)
ps4 = vpa(sol(19,1),3)
ps5 = vpa(sol(20,1),3)
ps6 = vpa(sol(21,1),3)
ps7 = vpa(sol(22,1),3)
ps8 = vpa(sol(23,1),3)
delt1 = vpa(sol(24,1),3)
delt2 = vpa(sol(25,1),3)
delt3 = vpa(sol(26,1),3)
delt4 = vpa(sol(27,1),3)
delt5 = vpa(sol(28,1),3)
delt6 = vpa(sol(29,1),3)
delt7 = vpa(sol(30,1),3)
p = vpa(sol(31,1),3)

% ------------------- dosazení vypočítaných hodnot --------------------

sigmaR1 = subs(sigmaR1)
sigmaR2 = subs(sigmaR2)
sigmaR3 = subs(sigmaR3)
sigmaR4 = subs(sigmaR4)
sigmaR5 = subs(sigmaR5)
sigmaR6 = subs(sigmaR6)
sigmaR7 = subs(sigmaR7)
sigmaR8 = subs(sigmaR8)

sigmaT1 = subs(sigmaT1)
sigmaT2 = subs(sigmaT2)
sigmaT3 = subs(sigmaT3)
sigmaT4 = subs(sigmaT4)
sigmaT5 = subs(sigmaT5)
sigmaT6 = subs(sigmaT6)
sigmaT7 = subs(sigmaT7)
sigmaT8 = subs(sigmaT8)

% sigmaZ1 = subs(sigmaZ1)
% sigmaZ2 = subs(sigmaZ2)
% sigmaZ3 = subs(sigmaZ3)
% sigmaZ4 = subs(sigmaZ4)
% sigmaZ5 = subs(sigmaZ5)
% sigmaZ6 = subs(sigmaZ6)
% sigmaZ7 = subs(sigmaZ7)
% sigmaZ8 = subs(sigmaZ8)


fplot(sigmaR1,[r1 r2],'color','red')
hold on
fplot(sigmaT1,[r1 r2],'color','red')
%fplot(sigmaZ1,[r1 r2],'color','red')

fplot(sigmaR2,[r2 r3],'color','blue')
fplot(sigmaT2,[r2 r3],'color','blue')
%fplot(sigmaZ2,[r2 r3],'color','blue')

fplot(sigmaR3,[r3 r4],'color','green')
fplot(sigmaT3,[r3 r4],'color','green')
%fplot(sigmaZ3,[r3 r4],'color','green')

fplot(sigmaR4,[r4 r5],'color','black')
fplot(sigmaT4,[r4 r5],'color','black')
%fplot(sigmaZ4,[r4 r5],'color','black')

fplot(sigmaR5,[r5 r6],'color','red')
fplot(sigmaT5,[r5 r6],'color','red')
%fplot(sigmaZ5,[r5 r6],'color','red')

fplot(sigmaR6,[r6 r7],'color','blue')
fplot(sigmaT6,[r6 r7],'color','blue')
%fplot(sigmaZ6,[r6 r7],'color','blue')

fplot(sigmaR7,[r7 r8],'color','green')
fplot(sigmaT7,[r7 r8],'color','green')
%fplot(sigmaZ7,[r7 r8],'color','green')

fplot(sigmaR8,[r8 r9],'color','black')
fplot(sigmaT8,[r8 r9],'color','black')
%fplot(sigmaZ8,[r8 r9],'color','black')

legend('sigmaR','sigmaT')


