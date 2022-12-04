clear all
clc
close all
% bezpecnost slabe zakriveneho prutu ohybaneho ve 3D
syms RA F1x F1y R fi x E Iy a G S F2 bet
% protoze jdeme od volneho konce smeremk vetknuti tak neni treba delat
% statickou rovnovahu, protoze jedina reakce, ktera bude v rovnicich
% figurovat je ta ve vzabe A, ktera je zaroven neznama. Tim je automaticky
% zaruceno, ze vsechny reakcni sily v rovnicich vnitrnich ucinku musi byt
% funkci te nezname rekace, podle ktere budeme nasledne derivovat


%Vnitrni ucinky prvni usek:

N1=-F1x*cos(fi)+F1y*sin(fi)
Ty1=F2
Tz1=F1x*sin(fi)+F1y*cos(fi)
Moy1=F1x*R*(1-cos(fi))+F1y*R*sin(fi)
%zde je ramenem kolma vzdalenost mezi silou F2 a nositelkou momentu Moz
%prochazejici stredem krivosti 
Moz1=-F2*R*sin(fi)
Mk1=-F2*R*(1-cos(fi))

%druhy usek

N2=F1x*sin(fi)-F1y*cos(fi)-RA*sin(fi)
Ty2=F2
Tz2=F1x*cos(fi)-F1y*sin(fi)-RA*cos(fi)
Moy2=F1x*R*(1+sin(fi))+F1y*R*cos(fi)-RA*R*sin(fi)
Moz2=-F2*R*cos(fi)
Mk2=-F2*R*(1+sin(fi))

%derivace vnitrnich ucinku dle RA. vidime ze nenulove budou jen ve druhem
%useku

DN2=diff(N2,RA);
DTz2=diff(Tz2,RA);
DMoy2=diff(Moy2,RA);

%deformacni podminka

uA=0==1/(E*Iy)*int(Moy2*DMoy2*R,fi,[0 pi/2])+bet/(G*S)*int(Tz2*DTz2*R,fi,[0 pi/2])+1/(E*S)*int(N2*DN2*R,fi,[0 pi/2])

E=70000;
bet=1.89;
d=35;
mu=0.3
G=E/(2*(1+mu));
S=pi*d^2/4;
R=300;

F1=2500;
F2=800;
alfa=15*pi/180;
F1x=F1*cos(alfa)
F1y=F1*sin(alfa)


%V�po�et kvadratick�ho momentu
Iy=pi*d^4/64;
%aby bylo Iy vy��sleno mus�me dosadit v�echny dosud zn�m� hodnoty (v na�em
%p��pad� b a h
Iy=subs(Iy);
%stejn� tak chceme dosadit zn�m� hodnoty do deformacni podminky
uA=subs(uA);

%nyn� m��eme vy�e�it rovnici a najit velikost reakce RA.
RA=vpa(solve(uA,RA),5)


% dosadime vse do relevantnich vnitrnich ucinku. tj nezajimaji nas
% posouvajici sily ktere genereuji max napeti radove menzi nez ohybove a kroutici momenty
%a navic uproistred prutu, kde jsou napeti od ohybu a krutu nulove uprostred prurezu

Moy1=subs(Moy1)
Moy2=subs(Moy2)
N1=subs(N1)
N2=subs(N2)
Moz1=subs(Moz1)
Moz2=subs(Moz2)
Mk1=subs(Mk1)
Mk2=subs(Mk2)

%jednotlive ucinky si vykreslime abychom videli jejich extremy
fplot(Moy1,[0 pi/2],'color','red')
hold on
fplot(Moy2,[0 pi/2],'color','blue')
hold off


fplot(Moz1,[0 pi/2],'color','red')
hold on
fplot(Moz2,[0 pi/2],'color','blue')
hold off

fplot(N1,[0 pi/2],'color','red')
hold on
fplot(N2,[0 pi/2],'color','blue')
hold off

fplot(Mk1,[0 pi/2],'color','red')
hold on
fplot(Mk2,[0 pi/2],'color','blue')
hold off

%pro redukovane napeti budou treba moduly prurezu v ohybu a krutu
Wo=Iy/(d/2)
Wk=(2*Iy)/(d/2)

%spocitame redukovane napeti. Absolutni hodnoty jsou nezbytne,
%protoze hodnoty normalovych sil i momentu nabyvaji jak kladnych tak
%zapornych hodnot a mohl by dojit k jejich odcitani.
%Druhy problem je prostorovy ohyb. diky prurezu tvaru kruhu je hlavni osou kazda osa prochazejici tezistem
%a s vyhodou tak lze VEKTOROVE secist oba momenty a tim dostaneme prosty
%ohyb vuci ose definovane orientaci tohoto vysledneho vektoru momentu
Mo1=sqrt(Moy1^2+Moz1^2)
Mo2=sqrt(Moy2^2+Moz2^2)
sigmaRED1=sqrt((abs(Mo1)/Wo+abs(N1)/S)^2+4*(abs(Mk1)/Wk)^2)
sigmaRED2=sqrt((abs(Mo2)/Wo+abs(N2)/S)^2+4*(abs(Mk2)/Wk)^2)
% vykreslime si prubehy redukovanych nap2ti
fplot(sigmaRED1,[0 pi/2],'color','red')
hold on
fplot(sigmaRED2,[0 pi/2],'color','blue')
hold off

%obrazek prubehu napeti dokazuje ze extrem napeti je na rozhrani intervalu
sigmaMAX=vpa(subs(sigmaRED1,fi,pi/2),5)
Re=280
kk=Re/sigmaMAX

% jak je videt, tak bezpecnost tohoto prutu je 1.2
%overeno ansysem. Vysledek v ansysu 243MPa, tedy chyba cca 5%
