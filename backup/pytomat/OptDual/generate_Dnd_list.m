function list=generate_Dnd_list(n,N_max,N_min,flag)
% list=generate_Dnd_list(n,N_max,N_min,flag)
% this function generates listing of torus indices with n unit cells under
% N_max*n number of atoms, 
% the output is a N-by-7 matrix, where the seven columns correspond to 
% (Natoms,m,g1,gm,z,E7T,E5T)
%
% if flag='off', then the stability criteria of TCNT is not enforced. This
% option is particularly useful in the list generation of buckled CNT.
%
% there are 9 kinds of even/odd sets possible combination of 4-indices
% (m,g1,gm,z), corresponding to total 16 kinds of (m,g1,gm,z,E7T,E5T) sets.
% see the following table:
% type | m g1 gm  z  |   E7T  E5T   |   minimal value
% ==========================================================               
%   1  | e  e  e  e  |   off  off   |  (2,0,0,2,0,0)                       
%      |             |   off   on   |  (2,0,2,2,0,1)                       
%      |             |    on  off   |  (2,2,0,2,1,0)                       
%      |             |    on   on   |  (2,2,2,2,1,1)                       
% ----------------------------------------------------------               
%   2  | o  o  e  o  |   off  off   |  (1,1,0,1,0,0)                       
%      |             |   off   on   |  (1,1,2,1,0,1)                       
% ----------------------------------------------------------               
%   3  | o  e  o  e  |   off  off   |  (1,0,1,2,0,0)                       
%      |             |    on  off   |  (1,2,1,2,1,0)                       
% ----------------------------------------------------------               
%   4  | e  e  o  e  |   off   on   |  (2,0,1,2,0,1)                       
%      |             |    on   on   |  (2,2,1,2,1,1)                       
% ----------------------------------------------------------               
%   5  | e  o  e  e  |    on  off   |  (2,1,0,2,1,0)                       
%      |             |    on   on   |  (2,1,2,2,1,1)                       
% ----------------------------------------------------------               
%   6  | e  o  o  o  |   off  off   |  (2,1,1,1,0,0)                       
% ----------------------------------------------------------               
%   7  | o  o  o  o  |   off   on   |  (1,1,1,1,0,1)                       
% ----------------------------------------------------------               
%   8  | o  o  o  e  |    on  off   |  (1,1,1,2,1,0)                       
% ----------------------------------------------------------               
%   9  | e  o  o  e  |    on   on   |  (2,1,1,2,1,1)                       
% ==========================================================               

if nargin<3
    N_min=0;
end
if nargin<4
    flag='';
end
list=zeros(0,7);

% case 1
temp1=pc_list(2,0,0,2,n,N_max);temp1=add_ET(temp1,0,0);
temp2=pc_list(2,0,2,2,n,N_max);temp2=add_ET(temp2,0,1);
temp3=pc_list(2,2,0,2,n,N_max);temp3=add_ET(temp3,1,0);
temp4=pc_list(2,2,2,2,n,N_max);temp4=add_ET(temp4,1,1);
list=[list;temp1;temp2;temp3;temp4]; clear temp1 temp2 temp3 temp4

% case 2
temp1=pc_list(1,1,2,1,n,N_max);temp1=add_ET(temp1,0,0);
temp2=pc_list(1,1,2,1,n,N_max);temp2=add_ET(temp2,0,1);
list=[list;temp1;temp2]; clear temp1 temp2

% case 3
temp1=pc_list(1,0,1,2,n,N_max);temp1=add_ET(temp1,0,0);
temp2=pc_list(1,2,1,2,n,N_max);temp2=add_ET(temp2,1,0);
list=[list;temp1;temp2]; clear temp1 temp2

% case 4
temp1=pc_list(2,0,1,2,n,N_max);temp1=add_ET(temp1,0,1);
temp2=pc_list(2,2,1,2,n,N_max);temp2=add_ET(temp2,1,1);
list=[list;temp1;temp2]; clear temp1 temp2

% case 5
temp1=pc_list(2,1,0,2,n,N_max);temp1=add_ET(temp1,1,0);
temp2=pc_list(2,1,2,2,n,N_max);temp2=add_ET(temp2,1,1);
list=[list;temp1;temp2]; clear temp1 temp2

% case 6
temp=pc_list(2,1,1,1,n,N_max);temp=add_ET(temp,0,0);
list=[list;temp]; clear temp

% case 7
temp=pc_list(1,1,1,1,n,N_max);temp=add_ET(temp,0,1);
list=[list;temp]; clear temp

% case 8
temp=pc_list(1,1,1,2,n,N_max);temp=add_ET(temp,1,0);
list=[list;temp]; clear temp

% case 9
temp=pc_list(2,1,1,2,n,N_max);temp=add_ET(temp,1,1);
list=[list;temp]; clear temp

% stability criterion
if ~strcmp(flag,'off')
    list(abs(list(:,3)-list(:,4))>=2*list(:,2),:)=[];
end

list(list(:,1)<=N_min*n,:)=[];
list=sortrows(list);
end

function temp=pc_list(m_min,g1_min,gm_min,z_min,n,N_max)
% parity-conserved list generating function
% n is the number of unit cells
% N_max is the upper bound of the total number of atoms

% N_max=N_max/n;
temp=zeros(0,5);
j=0;
for m=m_min:2:sqrt(N_max+4)-2
    for z=z_min:2:(N_max-m^2)/(4*m)
        for gm=gm_min:2:(N_max-4*m*z-2*m^2)/(2*(m+z))
            for g1=g1_min:2:(N_max-4*m*z-2*m^2-2*gm*(m+z))/(2*z)
                N = 2*g1*z + 2*(2*m*z+m^2) + 2*gm*(m+z);
                if N<=N_max
                    j=j+1;
                    temp(j,:)=[N m g1 gm z];
                end
            end
        end
    end
end
if numel(temp),temp(:,1)=temp(:,1)*n;end
end

function list=add_ET(list,E7T,E5T)

if E7T
    list=[list ones(size(list,1),1)];
else
    list=[list zeros(size(list,1),1)];
end

if E5T
    list=[list ones(size(list,1),1)];
else
    list=[list zeros(size(list,1),1)];
end
end
