function [R,g]=Dnd_TCNT_Dual(ind)

if nargin<1
    ind=[2 2 2 2 0 0];
end

m=ind(1);
g1=ind(2);
gm=ind(3);
z=ind(4);
g=m*2+g1+gm;
if ind(5)
    E7T=z/2;
    is=0;
else
    E7T=0;
    is=z/2;
end
if ind(6)
    E5T=(z+m)/2;
    os=0;
else
    E5T=0;
    os=(z+m)/2;
end
outer_SWT=0;
Dunlap_add=0;
R=Dn_tori(m,g1,gm,z,E7T,E5T,is,os,outer_SWT,Dunlap_add);
end
function [R,R_]=Dn_tori(m,g1,gm,z,E7T,E5T,inner_shift,outer_shift,outer_SWT,Dunlap_add)
% R=Dn_tori(m,g1,gm,z,E7T,E5T,inner_shift,outer_shift,outer_SWT,Dunlap_add)
% m is the distance between inner-ring 7-gons and outer-ring 5-gons
% g1 is the distance between inner-ring 7-gons
% gm is the distance between outer-ring 5-gons
% z is the length of the unit cell
% E7T is a transformation on the inner-ring, it can only take a value
%     between 1 and floor(z/2)
% E5T is a transformation on the outer-ring, it can only take a value
%     between 1 and floor((z+m)/2)
% inner_shift and outer_shift is a pair of numbers specifying the phase 
%     difference between inner-ring and outer-ring, Dnd and Dnh symmetry is
%     two extrema of the shifting parameters, it can only take a value 
%     between 0 and z-1. When the edge transformation parameter (E7T and 
%     E5T) is given, shifting parameters (inner_shift and outer_shift) will
%     no longer be free variables, else there would be 5-gons in the inner
%     ring or 7_gons in the outer ring.
% Given a set of (m,g1,gm,z), the number of atoms (N) is
%       N = 2*g1*z + 2*(2*m*z+m^2) + 2*gm*(m+z);

if nargin==0
    m=2;g1=2;gm=2;z=2;
end
if nargin<5
    E7T=0;
end
if nargin<6
    E5T=0;
end
if nargin<7
    inner_shift=-mod(g1,2)/2;
end
if nargin<8
    outer_shift=0;
end
if nargin<9
    outer_SWT=0;
end
if nargin<10
    Dunlap_add=0;
end

if outer_SWT&&(E5T||mod(gm,2)||(gm==0&&~Dunlap_add))
    disp 'outer_SWT is not allowed when E5T is on or gm is odd (gm>0)'
    disp 'outer_SWT is set to 0'
    outer_SWT=0;
end

if Dunlap_add&&((g1~=0||(gm~=0&&~outer_SWT)))
    disp 'Dunlap-type tori is allowd only when g1=0 and gm=0'
    disp 'g1 and gm are set to 0'
    g1=0;gm=0;
end

R_g1=g1_stripe(g1,z,E7T);
if m==0
    R_mid=[];
    R_mid_lower=[];
else
    R_mid=mid_stripe(z,m);
    if ~Dunlap_add
        R_mid_lower=[R_mid(:,1) -R_mid(:,2)+g1+2*m+gm];
        R_mid(:,2)=R_mid(:,2)+g1;
    else
        R_mid_lower=[R_mid(:,1) -R_mid(:,2)];
    end
end

if ~Dunlap_add
    R_gm=gm_stripe(gm,m,z,E5T,outer_SWT);
else
    R_gm=[];
end
if numel(R_gm)
    R_gm(:,2)=R_gm(:,2)+g1+m;
end


R=[R_g1;R_mid;R_mid_lower;R_gm];clear R_g1 R_mid R_mid_lower R_gm
if Dunlap_add
    R_add=Dunlap_addition(Dunlap_add,m,z,outer_SWT,gm);
    R=[R;R_add];
end

if E7T
    R(:,1)=shift(R,2*m+g1+gm,g1,g1/2);
    R(:,1)=shift(R,2*m+g1+gm,0,E7T);
    if inner_shift
        disp('inner shift is not a free variable when given a nonzero E7T value!')
    end
else
    if ~ispinteger(abs(inner_shift+g1/2),0)
        inner_shift=round(inner_shift)+mod(g1,2)/2;
        disp('2*inner_shift+g1 must be even')
        disp(['inner_shift has been setted to ' num2str(inner_shift)])
    end
    R(:,1)=shift(R,2*m+g1+gm,g1,inner_shift);
end

if E5T
    R(:,1)=shift(R,2*m+g1+gm,g1+m+gm,(-E5T-gm/2)*z/(z+m));
    if outer_shift
        disp('outer shift is not a free variable when given a nonzero E5T value!')
    end
else
    if ~ispinteger(abs(outer_shift+gm/2),0)
        outer_shift=round(outer_shift)+mod(gm,2)/2;
        disp('2*outer_shift+gm must be even')
        disp(['outer_shift has been setted to ' num2str(outer_shift)])
    end

    R(:,1)=shift(R,2*m+g1+gm,g1+m+gm,-outer_shift*z/(z+m));
end


if nargout==0
    dplot(R)
    set(gcf,'name',num2str([m g1 gm z E7T E5T inner_shift outer_shift outer_SWT Dunlap_add]),'numbertitle','off')
end

R(:,2)=R(:,2)-g1/2;

% RR=R;
% tt1=reshape(R(:,1),3,[]);
% ind=all(tt1<0,1);
% while any(ind)
%     tt1(:,ind)=tt1(:,ind)+z;
%     ind=all(tt1<0,1);
% end
% ind=all(tt1>z+Dunlap_add,1);
% while any(ind)
%     tt1(:,ind)=tt1(:,ind)-z-Dunlap_add;
%     ind=all(tt1>z+Dunlap_add,1);
% end
% R(:,1)=tt1(:);clear tt1 ind

R=1e-8*round(R*1e8);
R_=2*(z+Dunlap_add)+m;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=g1_stripe(g1,z,trans)
% R=g1_stripe(g1,z,trans)
% this function generate tiles in the first stripe (inner ring of the torus)
% trans=0 -> no transformation
% trans~=0 -> E7T, trans is a number between 1 and floor(z/2)
if ~g1
    R=[];
    return
end
if nargin<3
    trans=0;
end
R=[0 0;1 0;0.5 1;0 0;0.5 1;-0.5 1];
R=trans_expansion(0:g1-1,R,2);
temp=repmat(0:g1-1,6,1);
R(:,1)=R(:,1)-0.5*temp(:);clear temp
if ~trans
    R=trans_expansion(0:z-1,R,1);
elseif trans
    R=trans_expansion(0:z-trans-1,R,1);
    RR=[0 0;1 0;-0.5 1;1 0;0.5 1;-0.5 1];
    RR=trans_expansion(0:g1-1,RR,2);
    temp=repmat(0:g1-1,6,1);
    RR(:,1)=RR(:,1)-0.5*temp(:);clear temp
    RR=trans_expansion(0:trans-1,RR,1);RR(:,1)=RR(:,1)+z-trans;
    R=[R;RR];clear RR
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=mid_stripe(z,m)

temp0=[0 0;2 0;1 1;2 0;3 1;1 1];
temp0=trans_expansion(2*(0:z-1),temp0,1);
temp0=trans_expansion(0:m-1,temp0,2);
tt=repmat(0:m-1,6*z,1);temp0(:,1)=temp0(:,1)+tt(:);clear tt

temp1=[0 0;1 1;-1 1];
temp1=repmat(temp1,(m+1)*m/2,1);
temp1=temp1+tri_expansion(m-1)*[1 1;-1 1];
if m>1
    temp2=[-1 1;1 1;0 2];
    temp2=repmat(temp2,m*(m-1)/2,1);
    temp2=temp2+tri_expansion(m-2)*[1 1;-1 1];
else
    temp2=[];
end
R=[temp0;temp1;temp2];clear temp0 temp1 temp2

for k=1:m
    R(R(:,2)==k,1)=R(R(:,2)==k,1)+k;
    for l=2:2:2*(z+k)
        R((R(:,2)==k&R(:,1)==l),1)=R((R(:,2)==k&R(:,1)==l),1)-k*l/(z+k);
    end
end
R(:,1)=R(:,1)/2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=gm_stripe(gm,m,z,trans,outer_SWT)
if gm==0
    R=[];
    return
end
R=[0 0;1 0;0.5 1;1 0;0.5 1;1.5 1];
R(:,1)=R(:,1);
R=trans_expansion(0:gm-1,R,2);
temp=repmat(0:gm-1,6,1);
R(:,1)=R(:,1)+0.5*temp(:);clear temp
if ~trans
    R=trans_expansion((0:m+z-1),R,1);
    R=1e-3*round(R*1e3);
    if outer_SWT
        x=reshape(R(:,1),3,[]);
        y=reshape(R(:,2),3,[]);clear R
        cond1=all(y>=2*(z+m-x),1);
        cond2=all(y<=2*(x-z-m)+gm,1);
        cond3=all(y<=2*(z+m-x)+gm,1);
        x(:,cond1&cond2&cond3)=[];
        y(:,cond1&cond2&cond3)=[];
        R(:,1)=x(:);
        R(:,2)=y(:);
        temp=[0 0;-0.5 1;0 2;0 0;0.5 1;0 2];
        temp=trans_expansion((0:gm/2-1)/2,temp,1);
        tt=repmat(0:gm/2-1,6,1);temp(:,2)=temp(:,2)+tt(:);
        temp=trans_expansion(0:gm/2-1,temp,2);
        tt=-0.5*repmat(0:gm/2-1,3*gm,1);temp(:,1)=temp(:,1)+tt(:);
        R=[R;temp(:,1)+z+m temp(:,2)];
    end
else
    R=trans_expansion((0:m+z-trans-1),R,1);
    RR=[0 0;1 0;1.5 1;0 0;0.5 1;1.5 1];
    RR(:,1)=RR(:,1);
    RR=trans_expansion(0:gm-1,RR,2);
    temp=repmat(0:gm-1,6,1);
    RR(:,1)=RR(:,1)+0.5*temp(:);clear temp
    RR=trans_expansion((0:trans-1),RR,1);%RR(:,1)=RR(:,1)+z-trans;
%     figure;dplot(RR)
    R(:,1)=R(:,1)+trans;
    R=[R;RR];clear RR
%     figure;dplot(R)
end
R(:,1)=R(:,1)*z/(z+m);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Rx=shift(R,N,l,s)
% this function shifts the lth level of R, which is N in height, by s.
Rx=R(:,1);
R=[R (1:length(R))'];

tempy=reshape(R(:,2),3,[]);
I=repmat(any(tempy==l,1),3,1);
I_up=repmat(any(tempy==l,1) & all(tempy>=l,1),3,1);
I_down=repmat(any(tempy==l,1) & all(tempy<=l,1),3,1);

I=~I(:);
RR=[R(R(:,2)>l&I,1) R(R(:,2)>l&I,2)-l R(R(:,2)>l&I,3)
    R(R(:,2)<l&I,1) R(R(:,2)<l&I,2)+N-l R(R(:,2)<l&I,3)];

RR=[RR;R(I_up,1) R(I_up,2)-l R(I_up,3);R(I_down,1) R(I_down,2)+N-l R(I_down,3)];
[Ry,I]=sort(RR(:,2));
tt=[0 find([diff(Ry)' 1])];
RR(:,1)=RR(:,1)-s/2;
for j=1:N+1
    RR(I(tt(j)+1:tt(j+1)),1)=RR(I(tt(j)+1:tt(j+1)),1)+s*(j-1)/N;
end
Rx(RR(:,3))=RR(:,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R=Dunlap_addition(add,m,z,SWT,gm)

% R=[0 0;1 0;0 1;1 0;0 1;1 1;1 0;2 0;2 1;1 0;2 1;1 1];
R=[0 0;1 0.5;0 1;1 0.5;0 1;1 1.5;1 0.5;2 0;2 1;1 0.5;1 1.5;2 1];
R=trans_expansion((0:add-1)*2,R,1);
R=trans_expansion(0:2*m-1,R,2);

if SWT
    x=reshape(R(:,1),3,[]);
    y=reshape(R(:,2),3,[]);clear R
    cond1=all(y<=0.5*x,1);
    cond2=all(y<=-0.5*x+add,1);
    cond3=all(y>=-0.5*x+2*m,1);
    cond4=all(y>=0.5*x-add+2*m,1);
    cond=(cond1&cond2)|(cond3&cond4);
    x(:,cond)=[];
    y(:,cond)=[];
    R(:,1)=x(:);
    R(:,2)=y(:);
    
    temp=[0 0;1 -0.5;2 0;0 0;1 0.5;2 0];
    temp=trans_expansion(0:add-1,temp,1);
    tt=repmat(0:add-1,6,1);temp(:,2)=temp(:,2)+0.5*tt(:);
    temp=trans_expansion(0.5*(0:add-1),temp,2);
    tt=repmat(0:add-1,6*add,1);temp(:,1)=temp(:,1)-tt(:);
    temp(:,1)=temp(:,1)+add-1;
    temp(:,2)=temp(:,2)-0.5*(add-1);
    ind=all(reshape(temp(:,2),3,[])<=0,1);
    ind=repmat(ind,3,1);
    temp(ind(:),2)=temp(ind(:),2)+2*m;
    R=[R;temp];
end
if mod(gm,2)
    disp 'in this case, gm must be even'
    disp 'gm is set to 2'
    gm=2;
end
temp=[0 0;0 1;1 1;0 0;1 0;1 1;0 1;1 1;0 2;1 1;0 2;1 2];
temp1=trans_expansion(0:add-1,temp,1);
temp1=trans_expansion(2*(0:gm/2-1),temp1,2);
temp1(:,1)=temp1(:,1)*2;
temp1(:,2)=temp1(:,2)+2*m;
temp2=trans_expansion(0:z+m-1,temp,1);
temp2(:,1)=temp2(:,1)-z-m;
temp2=trans_expansion(2*(0:gm/2-1),temp2,2);
temp2(:,2)=temp2(:,2)+2*m;
R=[R;temp1;temp2];

R(:,1)=R(:,1)*z/(z+m)+z;
R(:,2)=R(:,2)-m;
end
