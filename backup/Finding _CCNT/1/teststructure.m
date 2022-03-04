function y=teststructure(x)

fileID3 = fopen('good.txt','r');
formatSpec2 = '%f';
ABCD = fscanf(fileID3,formatSpec2);

fileID35 = fopen('alreadydone.txt','r');
ABCD3 = fscanf(fileID35,formatSpec2);

fileID4 = fopen('bad.txt','r');
ABCD2 = fscanf(fileID4,formatSpec2);

ksd = sprintf('%1d',x);
iwant=str2double(ksd);

if ismember(iwant,ABCD) || ismember(iwant,ABCD2) || ismember(iwant,ABCD3)
    y=1;
else

try
   
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


     
     Min=min(pdist(XX)); %computing the minimum distance of all carbon atoms

     if Min<1.37         %condition that returns 1 (bad) when the distance between 2 C atoms are less than 1.37
        y=1;
     elseif isnan(radius)
        y=1;
     elseif isnan(XX(1,1))
        y=1;
     elseif x(1)==x(2)==x(3)==x(4)==1
        y=1;
     else
         y=0;
     end

catch
    y=1;
end

end
fclose all;
end


