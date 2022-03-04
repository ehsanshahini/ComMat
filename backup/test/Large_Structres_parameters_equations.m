% This file contains the geometrical paramteres and calculation of
% analytical stress vs MD vs modified.
% The first part solve the equations for stress and strain for pareto and
% give the pareto as a function of geometrical parameters

clear all;
close all;
clc;

syms zetta k alfa R d eps l sig;
assume(zetta>0 & k>0 & alfa>0 & R>0 & d>0 & l>0 & sig>0);

eq1=(sig==((d * l * eps * zetta * k)/(4*pi*R^2)));
eq2=(sig==7.449*eps^(-0.38));

S=solve(eq1,eq2,sig,eps);

sigmaY=S.sig;   %sigmaY as a function of geometrical parameteres
epsY=S.eps;     %%epsilonY as a function of geometrical parameteres



matsig=[];  % the matrice which contains stress from MD, Analytic, Modified, pitch angle

% %#1    This structure is way out of trend
% R = 29.730;
% N = 1.125;
% l0 = 99.3388;
% alfa = (12 * pi)/180;
% eps = 0.018;
% G = 43.23;
% EE = 115;
% d = 18.68;
% tanesh = 15.63;
% zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
% sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
% matsig=[matsig;tanesh,sig,sigmodified,alfa];

%#2
R = 29.128;
N = 1.1;
l0 = 73.80;
alfa = (9 * pi)/180;
eps = 2.691;
G = 43.23;
EE = 115;
d = 16.59;
SimulationNumber = sym(535);
tanesh = 3.527;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;  %The stress after modification
l=l0/N;
sigYfitted=eval(sigmaY); %the Stress from analytical equation independent from epsilon
epsYfitted=eval(epsY);   %the strain from analytical equation independent from epsilon
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#3
R = 25.55;
N = 1.165;
l0 = 116.38;
alfa = (16 * pi)/180;
eps = 1.56;
G = 43.23;
EE = 115;
d = 14.7285;
SimulationNumber = 525;
tanesh = 8.53;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#4
R = 24.795;
N = 1.25;
l0 = 184.25;
alfa = (25 * pi)/180;
eps = 0.862;
G = 43.23;
EE = 115;
d = 15.82;
tanesh = 9.687;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#5
R = 23.146;
N = 1.25;
l0 = 136.78;
alfa = (20 * pi)/180;
eps = 1.17;
G = 43.23;
EE = 115;
d = 15.186;
SimulationNumber = sym(89);
tanesh = 8.82;    
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#6
R = 17.149;
N = 1.5;
l0 = 274.38;
alfa = (40 * pi)/180;
eps = 0.288;
G = 43.23;
EE = 115;
d = 18.20;
tanesh = 10.18;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#7
index={[7,5,9,9,2,3]};
eps = 1.79;
tanesh = 6.35;
R = 28.563;
l0 = 118.69;
N = 1.25;
d = 18.85;
alfa = (13 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#8
tanesh = 6.970;
eps = 1.61;
index={[7,5,8,8,2,4]};
R = 25.692;
l0 = 115.42;
N = 1.25;
d = 17.52;
alfa = (14 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#9
tanesh = 5.54;
eps = 2.02;
index={[6,5,9,7,1,2]};
R = 24.895;
l0 = 78.49;
N = 1.25;
d = 16.13;
alfa = (12 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#10
tanesh = 12.18;
eps = 0.27;
index={[5,5,7,9,2,7]};
R = 15.606;
l0 = 282.36;
N = 1.5;
d = 16.64;
alfa = (44 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#11
tanesh = 13.16;
eps = 0.25;
index={[5,5,8,9,2,8]};
R = 15.65;
l0 = 290.35;
N = 1.5;
d = 17.55;
alfa = (46 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#12
tanesh = 9.59;
eps = 1.07;
index={[6,5,6,9,2,4]};
R = 24.413;
l0 = 164.81;
N = 1.25;
d = 15.18;
alfa = (23*pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#13
tanesh = 10.63;
eps = 0.28;
index={[6,5,8,9,2,9]};
R = 17.2630;
l0 = 281.75;
N = 1.5;
d = 19.98;
alfa = ( 43*pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];


%#14
tanesh = 15.38;
eps = 0.15;
index={[5,5,5,9,2,9]};
R = 10.584;
l0 = 312.66;
N = 1.75;
d = 17.535;
alfa = (56 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#15
tanesh = 3.65;
eps = 2.33;
index={[5,5,9,8,1,2]};
R = 28.531;
l0 = 75.42;
N = 1.0;
d = 15.90;
alfa = (11 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#16
tanesh = 3.60;
eps = 2.59;
index={[5,5,8,9,1,1]};
R = 30.239;
l0 = 75.31;
N = 1.1;
d = 16.00;
alfa = ( 9* pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#17
tanesh = 4.34;
eps = 2.12;
index={[7,5,9,8,1,3]};
R = 27.486;
l0 = 85.50;
N = 1.25;
d = 17.6285;
alfa = ( 11* pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#18
tanesh = 14.11;
eps = 0.21;
index={[5,5,6,8,2,9]};
R = 12.841;
l0 = 265.54;
N = 1.5;
d = 15.9532;
alfa = ( 46* pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#19
tanesh = 4.85;
eps = 2.085;
index={[5,5,9,7,1,1]};
R = 25.641;
l0 = 73.58;
N = 1.05;
d = 15.1721;
alfa = ( 12* pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#20
tanesh = 3.65;
eps = 2.19;
index={[6,5,9,8,2,2]};
R = 27.116;
l0 = 79.88;
N = 1.125;
d = 17.15;
alfa = (10 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#21
tanesh = 14.58;
eps = 0.20;
index={[5,5,7,8,2,9]};
R = 13.008;
l0 = 270.10;
N = 1.5;
d = 18.00;
alfa = (48 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#22
tanesh = 9.64;
eps = 0.91;
index={[5,5,5,9,2,3]};
R = 22.921;
l0 = 177.38;
N = 1.25;
d = 14.45;
alfa = ( 26* pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#23
tanesh = 10.33;
eps = 0.28;
index={[7,5,5,9,2,9]};
R = 16.860;
l0 = 272.74;
N = 1.5;
d = 19.37;
alfa = (40 * pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#24
tanesh = 13.78;
eps = 0.18;
index={[5,5,7,9,2,9]};
R = 12.431;
l0 = 308.74;
N = 1.5;
d = 17.78;
alfa = ( 50* pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];

%#25
tanesh = 2.98;
eps = 2.52;
index={[5,6,8,9,2,1]};
R = 29.861;
l0 = 72.81;
N = 1.05;
d = 14.89;
alfa = ( 9* pi)/180;
G = 43.23;
EE = 115;
zetta = 1/(cos(alfa)/(G * (1 + sin(alfa))) + (2 * sin(alfa) * tan(alfa))/(EE * (1 + sin(alfa))));
sig = (d * l0 * eps * zetta)/(4 * N * pi * R^2);
k=0.4442*exp(0.6661*alfa)-1.685*exp(-11.97*alfa);
sigmodified=sig*k;
l=l0/N;
sigYfitted=eval(sigmaY);
epsYfitted=eval(epsY);
matsig=[matsig;tanesh,sig,sigmodified,epsYfitted,sigYfitted,alfa];