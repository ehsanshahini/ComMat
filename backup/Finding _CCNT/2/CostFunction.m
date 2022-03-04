function [z,ytemp]=CostFunction(x,InputParameters)

global NFE;
if isempty(NFE)
    NFE=0;
end
NFE=NFE+1;

%% Helix
% function [XX,AA,t,u,radius,pt]=helix(type,index,shift,ncell)

index=[x(1) x(2) x(3) x(4) 0 0];
type=x(5);
shift=x(6);
ncell=15; % number of unit cells of the output


m=index(1);
g1=index(2);
gm=index(3);
z=index(4);
g=m*2+g1+gm;
if type==1 % Dnd unit cell
    if index(5)
        E7T=index(5);
        is=0;
    else
        E7T=0;
        is=z/2;
    end
    if index(6)
        E5T=index(6);
        os=0;
    else
        E5T=0;
        os=(z+m)/2;
    end
    outer_SWT=0;
    Dunlap_add=0;
elseif type==2 % Dnh unit cell
    E7T=0;
    E5T=0;
    os=0;
    is=0;
    outer_SWT=index(5);
    Dunlap_add=index(6);
end
R=Dn_tori(m,g1,gm,z,E7T,E5T,is,os,outer_SWT,Dunlap_add);
N=size(R,1)/3;

y=mean(reshape(R(:,2),3,[]));
R(reshape(repmat(y>m+gm+g1/2,3,1),[],1),2)=R(reshape(repmat(y>m+gm+g1/2,3,1),[],1),2)-g;
R(:,2)=R(:,2)+g1/2+m;
R(:,1)=R(:,1)-min(R(abs(R(:,2))<eps,1));
R(:,1)=R(:,1)*(z+m)/z;

r=norm([shift g])/2/pi;
th=atan2(shift,g);
u0=(z+m)*cos(th);
t0=(z+m)*sin(th)/r;

Rot=[cos(th) sin(th);-sin(th) cos(th)];
R=R*Rot;
Rcyl=zeros(size(R,1),3);
Rcyl(:,3)=R(:,1);
Rcyl(:,2)=R(:,2)/r;
Rcyl(:,1)=r;
R=round(R/Rot*1e5)/1e5;

Rxyz(:,1)=(Rcyl(:,1)+cos((R(:,2)-g+gm/2)/g*2*pi)/g*2*r).*cos(Rcyl(:,2));
Rxyz(:,2)=(Rcyl(:,1)+cos((R(:,2)-g+gm/2)/g*2*pi)/g*2*r).*sin(Rcyl(:,2));
Rxyz(:,3)=Rcyl(:,3);


Np=ceil(abs(shift)/(z+m))+1;
Rxyz=repmat(Rxyz,Np*2+1,1);
for k=1:Np
    Rxyz((k*2-1)*3*N+1:k*2*3*N,:)=bsxfun(@plus,Rxyz((k*2-1)*3*N+1:k*2*3*N,:),[0 0 -u0*k])*RotatM([t0*k 0 0 1]);
    Rxyz(k*2*3*N+1:(k*2+1)*3*N,:)=bsxfun(@plus,Rxyz(k*2*3*N+1:(k*2+1)*3*N,:),[0 0 u0*k])*RotatM([t0*k 0 0 -1]);
end
n=size(Rxyz,1)/N/3;

R=[
   bsxfun(@plus,R,[-(z+m) 0]);
   R;
   bsxfun(@plus,R,[(z+m) 0]);
  ];
for k=2:Np
    R=[R;
       bsxfun(@plus,R(N*3+1:N*3*2,:),[-(z+m)*k 0]);
       bsxfun(@plus,R(N*3+1:N*3*2,:),[(z+m)*k 0]);
      ];
end

% making the shift
R(:,1)=R(:,1)-R(:,2)/g*shift;

A=conn_c(round(Rxyz*1e5));

Rxyz=squeeze(mean(reshape(Rxyz,3,[],3),1));

[m1,m2]=find(A);
D=mean(sqrt(sum((Rxyz(m1,:)-Rxyz(m2,:)).^2,2)));
Rxyz=Rxyz/D*1.43; % averaging the bond lengths
u0=u0/D*1.43;
clear m1 m2 D

A=conn_c([mod(round(R(:,1)*1e5)/1e5,n*(z+m)) mod(round(R(:,2)*1e5)/1e5,g)]);
Rxyz=Rxyz+rand(size(Rxyz))*1e-3;

% optimization routine
[X,t,u]=BFGS_helix(Rxyz,A,n,t0,u0,'');


% outputting: generating the desired length of the helix
XX=[
    bsxfun(@plus,X,[0 0 -u])*RotatM([t 0 0 1]);
    X;
    bsxfun(@plus,X,[0 0 u])*RotatM([t 0 0 -1]);
    ];
for k=2:ncell/2
    XX=[XX;
        bsxfun(@plus,XX(N+1:N*2,:),[0 0 -u*k])*RotatM([t*k 0 0 1]);
        bsxfun(@plus,XX(N+1:N*2,:),[0 0 u*k])*RotatM([t*k 0 0 -1]);
        ];
end
R=[
   bsxfun(@plus,R(1:N*3,:),[-(z+m) 0]);
   R(1:N*3,:);
   bsxfun(@plus,R(1:N*3,:),[(z+m) 0]);
  ];
for k=2:ncell/2
    R=[R;
       bsxfun(@plus,R(N*3+1:N*3*2,:),[-(z+m)*k 0]);
       bsxfun(@plus,R(N*3+1:N*3*2,:),[(z+m)*k 0]);
      ];
end
AA=conn_c([R(:,1) mod(round(R(:,2)*1e5)/1e5,g)]);

% calculate structure parameters:
radius=mean([sqrt(max(X(:,1).^2+X(:,2).^2)),sqrt(min(X(:,1).^2+X(:,2).^2))]);
pt=atan2(abs(u),radius*t)/pi*180;



%% Data File LAMMPS 

if isfile('data.cnt')
   delete('data.cnt');
end
Min = min(XX,[],1);
Max = max(XX,[],1);
MinDis=min(pdist(XX));

if MinDis>1.37
%calculating the length and arc length of the CCNT
Length=Max(3)-Min(3);
ArcL=3.2*sqrt((2*pi*radius)^2+((Length/3.2)^2));
adding=30+((ArcL-Length)/2);
fid = fopen('data.cnt', 'w');

fprintf(fid,'The data file data.cnt:\n\n');
fprintf(fid,'%d atoms\n',length(XX));
fprintf(fid,'1 atom types\n');
fprintf(fid, '%0.1f %0.1f      xlo xhi\n',Min(1)-30,Max(1)+30);
fprintf(fid, '%0.1f %0.1f      ylo yhi\n',Min(2)-30,Max(2)+30);
fprintf(fid, '%0.1f %0.1f      zlo zhi\n\n',Min(3)-adding,Max(3)+adding);
fprintf(fid, 'Atoms\n\n');

for i=1:length(XX)
    fprintf(fid,'%i	1	%f	%f	%f\n',i,XX(i,1),XX(i,2),XX(i,3));
end
fclose(fid);

%% Structure Parameters

diameter = sqrt(max(X(:,1).^2+X(:,2).^2))-sqrt(min(X(:,1).^2+X(:,2).^2));
lengthccnt = Max(1,3)-Min(1,3);

xcurve = XX(:,1);
zcurve = XX(:,3);

[xData, yData] = prepareCurveData( zcurve, xcurve );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

% Fit model to data.
[fitresult, ~] = fit( xData, yData, ft, opts );

syms zl;
assume(Min(3) < zl < Max(3));

a1=fitresult.a1;
b1=fitresult.b1;
c1=fitresult.c1;

try
    eq = a1*sin(b1*zl+c1);
    difeq = diff(eq);
    Solvation = sort(eval(solve(difeq,zl)));

    lengthSolv = length(Solvation)/2;
    if length(Solvation) < 2
    pitchl = lengthccnt;
    else
        indice1 = ceil(lengthSolv);
        indice2 = indice1+1;
        pitchl = 2*abs(Solvation(indice1)-Solvation(indice2));
    end

    numberofturns = lengthccnt / pitchl;
    pitchangle = 180/pi*(atan(pitchl/(2*pi*radius)));
catch
    pitchl = 0; numberofturns = 0; pitchangle = 0;
end

inradius = InputParameters(1); 
inpitchangle = InputParameters(2);
indiameter = InputParameters(3);


radiusC = abs(radius - inradius);
pitchangleC = abs(pitchangle - inpitchangle);
diameterC = abs(diameter - indiameter);


z=[radiusC
   pitchangleC
   diameterC];    %***Check***%

%% Put all the results into ytemp

ytemp.Position=x;
ytemp.SimulationNo=NFE;
ytemp.radius=radius;
ytemp.diameter=diameter;
ytemp.lengthccnt=lengthccnt;
ytemp.pitchl=pitchl;
ytemp.numberofturns=numberofturns;
ytemp.pitchangle=pitchangle;

else
    z=[0;0];
    ytemp.Position=x;
    ytemp.SimulationNo=NFE;
    ytemp.radius=0;
    ytemp.diameter=0;
    ytemp.lengthccnt=0;
    ytemp.pitchl=0;
    ytemp.numberofturns=0;
    ytemp.pitchangle=0;
    
end
fclose all;
end
