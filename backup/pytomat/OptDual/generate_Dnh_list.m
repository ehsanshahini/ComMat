function list=generate_Dnh_list(n,N_max,N_min,flag)
% list=generate_Dnh_list(n,N_max,N_min,flag)
% this function generates list of (stable) Dnh-symmetric TCNT indices set
% with rotational symmetry number n between N_max and N_min atoms per unit
% cell.
%
% if flag='off', then the stability criteria of TCNT is not enforced. This
% option is particularly useful in the list generation of buckled CNT.
%
% existence criteria for Dnh TCNT:
% if irT=0, then g1 must be even
% if irT=1, then z must be even
% if orT=0, then gm must be even
% if orT=1; then z+m must be even and the smallest possible value for gm is
%           (z-m)/2
% if orT=2; then gm must be a nonzero even number
% 1. if (irT,orT)=(0,0), then g1 and gm must be even
% 2. if (irT,orT)=(0,1), then g1 and z+m must be even and the smallest
%                        possible value for gm is (z-m)/2
% 3. if (irT,orT)=(1,0), then gm and z must be even 
% 4. if (irT,orT)=(1,1), then z and m must be even and the smallest
%                        possible value for gm is (z-m)/2
% 5. if (irT,orT)=(0,2), then both g1 and gm must be even, and gm cannot be
%                        zero; and gm must be smaller than 2*(z+m),
%                        otherwise it reduces to the case
%                        (2*m+z,g1,gm/2-z-m,z,0,1)
% 6. if (irT,orT)=(1,2), then z and gm must be even and gm cannot be zero

% ============ to be done ==================================
% there are 6 kinds of even/odd sets possible combination of 4-indices
% (m,g1,gm,z), corresponding to total 13 kinds of (m,g1,gm,z,irT,orT) sets.
% see the following table:
% type | m g1 gm  z  |   irT  orT   |   minimal value
% ==========================================================               
%   1  | e  e  e  e  |     0    0   |  (2,0,0,2,0,0)                       
%      |             |     0    1   |  (2,0,2,2,0,1)                       
%      |             |     1    0   |  (2,0,0,2,1,0)                       
%      |             |     1    1   |  (2,0,2,2,1,1)                       
%      |             |     0    2   |  (2,0,2,2,0,2)                       
% ----------------------------------------------------------               
%   2  | e  e  e  o  |     0    0   |  (2,0,0,1,0,0)                       
%      |             |     0    2   |  (2,0,2,1,0,2)                       
% ----------------------------------------------------------               
%   3  | e  e  o  e  |     0    1   |  (2,0,1,2,0,1)                       
%      |             |     1    1   |  (2,0,1,2,1,1)                       
% ----------------------------------------------------------               
%   4  | e  o  e  e  |     1    0   |  (2,1,0,2,1,0)                       
%      |             |     1    1   |  (2,1,0,2,1,1)                       
%      |             |     1    2   |  (2,1,2,2,1,2)                       
% ----------------------------------------------------------               
%   5  | e  o  o  e  |     1    1   |  (2,1,1,2,1,1)                       
% ----------------------------------------------------------               
%   6  | o  e  e  e  |     0    0   |  (3,0,0,2,0,0)                       
%      |             |              |  (1,2,0,2,0,0)                       
%      |             |              |  (1,0,2,2,0,0)                       
%      |             |     1    0   |  (3,0,0,2,1,0)                       
%      |             |     0    2   |  (1,0,2,2,0,2)                       
% ----------------------------------------------------------               
% ==========================================================               



if nargin<1
    n=5;
end
if nargin<2
    N_max=40;
end
if nargin<3
    N_min=0;
end
if nargin<4
    flag='';
end

list=[];

% 1. (irT,orT)=(0,0) case
list1=[];
j=0;
for m=1:sqrt(N_max+4)-2
    for z=1:(N_max-m^2)/(4*m)
        for gm=0:2:(N_max-4*m*z-2*m^2)/(2*(m+z))
            for g1=0:2:(N_max-4*m*z-2*m^2-2*gm*(m+z))/(2*z)
                N = 2*(2*m*z+m^2) + 2*g1*z + 2*gm*(m+z);
                if N<=N_max
                    j=j+1;
                    list1(j,:)=[N m g1 gm z];
                end
            end
        end
    end
end
list1=[list1 zeros(size(list1,1),2)];
if numel(list1),list=[list;list1];end


% 2. (irT,orT)=(0,1) case
list2=[];
j=0;
for m=1:sqrt(N_max*2/3)
    for z=1:sqrt(abs(12*m^2-2*N_max))+3*m
        for gm=0:(N_max-3/2*m^2-3*m*z+z^2/2)/(2*(m+z))
            for g1=0:2:(N_max-3/2*m^2-3*m*z+z^2/2-(4*gm-z-m)*(z+m)/2)/(2*z)
                N = 2*(2*m*z+m^2) + 2*g1*z + (m+z)/2*(4*gm-z-m);
                if N<=N_max&&~mod(z+m,2)&&N>0
                    j=j+1;
                    list2(j,:)=[N m g1 gm z];
                end
            end
        end
    end
end
% the smallest possible for gm is (z-m)/2
list2(list2(:,4)<(list2(:,5)-list2(:,2))/2,:)=[]; 
% girth too small
list2(list2(:,2)==2&list2(:,3)==0&list2(:,4)<list2(:,5)/2,:)=[];

list2=[list2 zeros(size(list2,1),1) ones(size(list2,1),1)];
if numel(list2),list=[list;list2];end

% 3. (irT,orT)=(1,0) case
list3=[];
j=0;
for m=1:sqrt((N_max+6)/2)-2
    for z=2:sqrt((N_max+6*m^2)*2)-m*4
        for gm=0:2:(N_max-2*m*(m+2*z)-z^2/2)/(2*(m+z))
            for g1=0:(N_max-2*m*(m+2*z)-z^2/2-2*gm*(z+m))/(2*z)
                N = 2*(2*m*z+m^2) + z/2*(z+4*g1) + 2*gm*(z+m);
                if N<=N_max&&~mod(z,2)
                    j=j+1;
                    list3(j,:)=[N m g1 gm z];
                end
            end
        end
    end
end
list3=[list3 ones(size(list3,1),1) zeros(size(list3,1),1)];
if numel(list3),list=[list;list3];end


% 4. (irT,orT)=(1,1) case
list4=[];
j=0;
for m=1:sqrt(N_max*2/3)
    for z=1:(N_max-3/2*m^2)/3/m
        for gm=0:(N_max-3/2*m^2-3*m*z)/(2*(m+z))
            for g1=0:(N_max-3/2*m^2-3*m*z-2*gm*(z+m))/(2*z)
                N = 2*(2*m*z+m^2) + z/2*(z+4*g1) + (z+m)/2*(4*gm-z-m);
                if N<=N_max&&~mod(z,2)&&~mod(z+m,2)
                    j=j+1;
                    list4(j,:)=[N m g1 gm z];
                end
            end
        end
    end
end
% the smallest possible for gm is (z-m)/2
list4(list4(:,4)<(list4(:,5)-list4(:,2))/2,:)=[];
% girth too small
list4(list4(:,2)==2&list4(:,3)==0&list4(:,4)==0,:)=[];
list4=[list4 ones(size(list4,1),2)];
if numel(list4),list=[list;list4];end

% 5. (irT,orT)=(0,2) case
list5=list1;
list5(~list5(:,4),:)=[]; % gm>0
list5(list5(:,4)/2>=(list5(:,2)+list5(:,5)),:)=[]; % gm<2*(z+m)
list5(:,6:7)=repmat([0 2],size(list5,1),1);
if numel(list5),list=[list;list5];end

% 6. (irT,orT)=(1,2) case
list8=list3;
list8(~list8(:,4),:)=[]; % gm>0
list8(list8(:,4)/2>=(list8(:,2)+list8(:,5)),:)=[]; % gm<2*(z+m)
list8(:,6:7)=repmat([1 2],size(list8,1),1);
if numel(list8),list=[list;list8];end

if numel(list),list(:,1)=list(:,1)*n;end

% stability criterion
if ~strcmp(flag,'off')
    list(list(:,3)./sum(list(:,[2 2 3 4]),2)>=0.5,:)=[];
    list(list(list(:,7)~=2,4)./sum(list(list(:,7)~=2,[2 2 3 4]),2)>=0.5,:)=[];
    list(sum(list(:,[2 2 3 4]),2)<3,:)=[];
    gmt=(list(:,4)-(list(:,5)+list(:,2))/2);gmt=gmt.*(gmt>0);
    list(list(:,2)*2+list(:,3)+gmt<3&list(:,7)==1,:)=[];
end

list(list(:,1)<=N_min*n,:)=[];
list=sortrows(list);
end