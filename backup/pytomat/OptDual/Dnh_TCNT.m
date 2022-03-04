function [R_xyz,A]=Dnh_TCNT(ind,n)
% [R_xyz,A]=Dnh_TCNT(ind,n)

% close all
% clear all
% ind=[2 12 10 1 0 0];
% n=7;

if nargin<1
    ind=[2 2 2 2 0 0];
end
if nargin<2
    n=5;
end

m=ind(1);g1=ind(2);gm=ind(3);z=ind(4);
E7T=ind(5);
switch ind(6)
    case 0
        E5T=0;
        oSWT=0;
    case 1
        E5T=1;
        oSWT=0;
    case 2
        E5T=0;
        oSWT=1;
end

if mod(g1,2)&&~E7T
    warning(['g1 and must be even when E7T is NOT enabled!' 10 'g1 is set to g1+1'])
    g1=g1+1;
end
if mod(gm,2)&&~E5T
    warning(['gm and must be even when E5T is NOT enabled!' 10 'gm is set to gm+1'])
    gm=gm+1;
end
if E7T&&mod(z,2)
    warning(['try an even z' 10 'z is set to z+1'])
    z=z+1;
end
if E5T&&mod(z+m,2)
    warning(['try an even (z+m)' 10 'm is set to m+1'])
    m=m+1;
end
if E5T&&gm<(z-m)/2
    warning(['when E5T is on, the smallest possible value for gm is (z-m)/2!' ...
             10 'gm is set to (z-m)/2 = ' num2str((z-m)/2)])
    gm=(z-m)/2;
end
if oSWT&&E5T
    warning(['oSWT is not allowed when E5T is on!' 10 'switch off oSWT'])
    oSWT=0;
end
if gm==0&&oSWT
    warning(['when gm=0, oSWT is not allowed' 10 'switch of oSWT'])
    oSWT=0;
end
if oSWT&&gm/2>=(z+m)
    warning(['when gm/2>=(z+m), oSWT reduces to the case (2*m+z,g1,gm/2-z-m,z,0,1,0,0)' 10 'switch off oSWT'])
    oSWT=0;
end

if E5T
    gmt=(gm-(z+m)/2)*((gm-(z+m)/2)>0);
else
    gmt=gm;
end

g=m*2+g1+gmt;
if g<3
    error(['the girth of the TCNT is too small!' 10 'try indices of larger values'])
elseif m==2&&gm==0&&E5T&&g1==0
    error(['the girth of the TCNT is too small!' 10 'try indices of larger values'])
end

Natoms=2*m*(m+2*z);
if E7T
    Natoms=Natoms+z/2*(z+4*g1);
else
    Natoms=Natoms+2*g1*z;
end
if E5T
    Natoms=Natoms+(z+m)/2*(4*gm-z-m);
else
    Natoms=Natoms+2*gm*(m+z);
end
% Natoms=Natoms*n;

Rm=mid_stripe(z,m);
if ~gm&&oSWT
    ind=ceil(find(abs(Rm(:,1))<eps&abs(Rm(:,2)-m)<eps)/3)*3;
    Rm(ind-2:ind,:)=[];
    Rm(abs(Rm(:,1)-z)<eps&abs(Rm(:,2)-m)<eps,2)=...
    Rm(abs(Rm(:,1)-z)<eps&abs(Rm(:,2)-m)<eps,2)+1;
end

if E7T
    if g1
        Rg1=mid_stripe(g1,z/2);
        Rg1=[z/2-Rg1(:,2) Rg1(:,1)];
        Rg1=[Rg1;z-Rg1(:,1) Rg1(:,2)];
    else
        Rg1=g0_tri(z/2);
        Rg1=[Rg1;z-Rg1(:,1) Rg1(:,2)];
        Rm(:,1)=Rm(:,1)-z/2;
        Rm(abs(Rm(:,2))<eps,2)=abs(Rm(abs(Rm(:,2))<eps,1)/z);
        Rm(:,1)=Rm(:,1)+z/2;
    end
else
    Rg1=g_stripe(g1,z);
end
if E5T
    if gm>(z+m)/2
        Rgm=mid_stripe(gm-(z+m)/2,(z+m)/2);
        Rgm=[Rgm(:,2)*z/(z+m) Rgm(:,1)];
        Rgm=[Rgm;z-Rgm(:,1) Rgm(:,2)];
    elseif gm==(z+m)/2
        Rgm=g0_tri((z+m)/2);
        Rgm(:,1)=z/2-Rgm(:,1)/(z+m)*z;
        Rgm=[Rgm;z-Rgm(:,1) Rgm(:,2)];
        Rm=bsxfun(@minus,Rm,[z/2 m]);
        Rm(abs(Rm(:,2))<eps,2)=abs(Rm(abs(Rm(:,2))<eps,1)/z)-0.5;
        Rm=bsxfun(@plus,Rm,[z/2 m]);
    else
        tt=(z+m)/2-gm;
        if tt>m
            gt=gm;
            gm=(z-m)/2;
        end
        k=(z+m)/2-gm;
        Rm(:,1)=Rm(:,1)-z/2;
        Rm(:,1)=Rm(:,1).*(Rm(:,2)+z)/z;
        x=mean(reshape(Rm(:,1),3,[]));
        y=mean(reshape(Rm(:,2),3,[]));
        ind=reshape(repmat(y>-2*x+2*gm+m,3,1),[],1);
        Rm(ind,:)=[];
        x=mean(reshape(Rm(:,1),3,[]));
        y=mean(reshape(Rm(:,2),3,[]));
        ind1=reshape(repmat(y>2*x+2*gm+m,3,1),[],1);
        ind2=abs(Rm(:,2)-m)<eps&~ind1;
        Rm(ind1,:)=bsxfun(@plus,Rm(ind1,:),[gm -m]);
        Rm(ind1,1)=Rm(ind1,1)-Rm(ind1,2)/2;
        Rm(ind1,:)=Rm(ind1,:)*[0.5 -1;0.5 1];
        Rm(ind1,:)=bsxfun(@minus,Rm(ind1,:),[gm -m]);
        Rm(Rm(:,2)<(m-z)/2+gm,1)=Rm(Rm(:,2)<(m-z)/2+gm,1)./(Rm(Rm(:,2)<(m-z)/2+gm,2)+z)*z;
        Rm(Rm(:,2)>=(m-z)/2+gm,1)=Rm(Rm(:,2)>=(m-z)/2+gm,1)/(m+z+gm*2)*z*2;
        if gm,Rm(ind2,2)=Rm(ind2,2)-0.5+abs(Rm(ind2,1))*(m+z+gm*2)/z/gm/4;end
        Rm(:,1)=Rm(:,1)+z/2;

        Rgm=g0_tri(gm);
        Rgm(:,1)=z/2-Rgm(:,1)/(m+z+gm*2)*z*2;
        Rgm=[Rgm;z-Rgm(:,1) Rgm(:,2)];
        if tt>m
            gm=gt;
        end
    end
else
    Rgm=g_stripe(gm,z+m);
    if oSWT&&gm
        x=mean(reshape(Rgm(:,1),3,[]));
        y=mean(reshape(Rgm(:,2),3,[]));
        ind=reshape(repmat(y>2*x&y<-2*x+gm,3,1),[],1);
        Rgm(ind,1)=Rgm(ind,1)+z+m;
        Rgm=bsxfun(@minus,Rgm,[z+m gm/2]);
        x=mean(reshape(Rgm(:,1),3,[]));
        y=mean(reshape(Rgm(:,2),3,[]));
        ind=reshape(repmat(y<2*x+gm/2&-2*x+gm/2&y>2*x-gm/2&y>-2*x-gm/2,3,1),[],1);
        Rgm(ind,:)=[Rgm(ind,2)/2 Rgm(ind,1)*2];
        Rgm=bsxfun(@plus,Rgm,[z+m gm/2]);
        if gm/2>(z+m)
            x=mean(reshape(Rgm(:,1),3,[]));
            y=mean(reshape(Rgm(:,2),3,[]));
            ind=reshape(repmat(y<2*x-(z+m)*4+gm&y>-2*x+2*(z+m),3,1),[],1);
            Rgm(ind,:)=bsxfun(@minus,Rgm(ind,:),[1 2]*(gm/2-z-m)/2);
        end
    end
    if gm,Rgm(:,1)=Rgm(:,1)*z/(z+m);end
end

if numel(Rgm),Rgm=[Rgm(:,1) Rgm(:,2)+g1+m];end

R=[Rg1;Rm(:,1) Rm(:,2)+g1;Rgm;z-Rm(:,1) -Rm(:,2)+g1+m*2+gmt];
R_xyz=TCNT_int_geom(R,[m g1 gmt z],n);
% A=conn_c(round(R_xyz*1e5));
A=conn(round(R_xyz*1e5));
[m1,m2]=find(A);
R_xyz=squeeze(mean(reshape(R_xyz,3,[],3)));
D=mean(sqrt(sum((R_xyz(m1,:)-R_xyz(m2,:)).^2,2)));
R_xyz=R_xyz/D*1.43;
end