function [R_xyz,A]=Dnd_TCNT(ind,n)
% [R_xyz,A]=Dnd_TCNT(ind,n)

% clear all
% close all
% ind=[1 2 1 1 0 1];
% n=5;

if nargin<1
    ind=[2 2 2 2 0 0];
end
if nargin<2
    n=5;
end

m=ind(1);g1=ind(2);gm=ind(3);z=ind(4);
if length(ind)<6
    ind=[ind 0 0];
end
E7T=ind(5);E5T=ind(6);is=z/2;os=(z+m)/2;
if ~mod(g1,2) && mod(z,2)
    error('When g1 is even and z is odd, Dnd does not exist.')
end
if mod(g1,2) && ~mod(z,2) && ~E7T
    disp('E7T must be enabled when g1 is odd and z is even')
    E7T=z/2;
end
if mod(g1,2) && mod(z,2) && E7T
    disp('E7T must be disabled when g1 is odd and z is odd')
    E7T=0;
end
if ~mod(gm,2) && mod(m+z,2)
    error('When gm is even and z+m is odd, Dnd does not exist.')
end
if mod(gm,2) && ~mod(z+m,2) && ~E5T
    disp('E5T must be enabled when gm is odd and z+m is even')
    E5T=(z+m)/2;
end
if mod(gm,2) && mod(z+m,2) && E5T
    disp('E5T must be disabled when gm is odd and z+m is odd')
    E5T=0;
end
if ~g1&&E7T
    disp 'E7T makes no difference when g1=0'
    E7T=0;
end
if ~gm&&E5T
    disp 'E5T makes no difference when gm=0'
    E5T=0;
end

if E7T
    E7T=z/2;
    is=0;
end
if E5T
    os=0;
    E5T=(z+m)/2;
end
R=Dn_tori(m,g1,gm,z,E7T,E5T,is,os);
R(:,2)=R(:,2)+g1/2;
R_xyz=TCNT_int_geom(R,[m g1 gm z],n);
A=conn(round(R_xyz*1e5));
% A=conn_c(round(R_xyz*1e5));
[m1,m2]=find(A);
R_xyz=squeeze(mean(reshape(R_xyz,3,[],3)));
D=mean(sqrt(sum((R_xyz(m1,:)-R_xyz(m2,:)).^2,2)));
R_xyz=R_xyz/D*1.43;
end