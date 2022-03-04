clear;
clc;
A=struct([]);
YoungMod=[];
for i=1:5
    ii=num2str(i);
    cd(ii);
    s=dlmread('stress_strain.txt');
    A(i).ss=s;
    cd ..;
    x=s(:,1);
    y=s(:,2);
    f=fit(x,y,'poly1');
    E=f.p1;
    YoungMod=[YoungMod;E];
 end
