function Parameters=helix(type,index,shift)
try
index=[index,0,0];
% type=1;
% shift=1;
ncell=15; % number of unit cells of the output


m=index(1);   %n75
g1=index(2);  %n77
gm=index(3);  %n55
z=index(4);   %s
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
th=atan2(shift,g);   % it has a relation with pitch angle
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

% figure(2);
% plot3(XX(:,1),XX(:,2),XX(:,3),'o','MarkerSize',5);

%% Data File LAMMPS 

Min = min(XX,[],1);
Max = max(XX,[],1);

fid = fopen('data.cnt', 'w');
fprintf(fid,'The data file data.cnt:\n\n %d atoms\n 1 atom types\n %0.1f %0.1f xlo xhi\n %0.1f %0.1f ylo yhi\n %0.1f %0.1f zlo zhi\n\nAtoms\n\n',length(XX),Min(1)-30,Max(1)+30,Min(2)-30,Max(2)+30,Min(3)-80,Max(3)+80);
for i=1:length(XX)
    fprintf(fid,'%i	1	%f	%f	%f\n',i,XX(i,1),XX(i,2),XX(i,3));
end

Min2 = min(X,[],1);
Max2 = max(X,[],1);

fid = fopen('dataX.cnt', 'w');
fprintf(fid,'The data file data.cnt:\n\n %d atoms\n 1 atom types\n %0.1f %0.1f xlo xhi\n %0.1f %0.1f ylo yhi\n %0.1f %0.1f zlo zhi\n\nAtoms\n\n',length(X),Min2(1),Max2(1),Min2(2),Max2(2),Min2(3),Max2(3));
for i=1:length(X)
    fprintf(fid,'%i	1	%f	%f	%f\n',i,X(i,1),X(i,2),X(i,3));
end

fclose('all');

%% Structure Parameters and curve fitting

radius=mean([sqrt(max(X(:,1).^2+X(:,2).^2)),sqrt(min(X(:,1).^2+X(:,2).^2))]);
pt=atan2(abs(u),radius*t)/pi*180;
diameter = sqrt(max(X(:,1).^2+X(:,2).^2))-sqrt(min(X(:,1).^2+X(:,2).^2));
lengthccnt = Max(1,3)-Min(1,3);

xcurve = XX(:,1);
ycurve = XX(:,2);
zcurve = XX(:,3);
try

[xData, yData] = prepareCurveData( zcurve, xcurve );

% Set up fittype and options.
ft = fittype( 'sin1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
 figure( 'Name', 'xz-plane' );
 h = plot( fitresult, xData, yData );

syms zl;
assume(Min(3) < zl < Max(3));

a1=fitresult.a1;
b1=fitresult.b1;
c1=fitresult.c1;

    eq = a1*sin(b1*zl+c1);
    difeq = diff(eq);
    Solvation = sort(eval(solve(difeq,zl)));

    lengthSolv = length(Solvation)/2;
    if length(Solvation) < 2
    Pitchl = lengthccnt;
    else
        indice1 = ceil(lengthSolv);
        indice2 = indice1+1;
        Pitchl = 2*abs(Solvation(indice1)-Solvation(indice2));
    end

    numberofturns = lengthccnt / Pitchl;
    pitchangle = 180/pi*(atan(Pitchl/(2*pi*radius)));
catch
    Pitchl = 0; numberofturns = 0; pitchangle = 0; radius = 0; diameter = 0; lengthccnt =0;
end
Parameters = [radius,diameter,pitchangle,Pitchl,numberofturns,lengthccnt];
catch
    Pitchl = 0; numberofturns = 0; pitchangle = 0; radius = 0; diameter = 0; lengthccnt =0;
    Parameters = [radius,diameter,pitchangle,Pitchl,numberofturns,lengthccnt];
end