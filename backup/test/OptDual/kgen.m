function [ii,k]=kgen(A,sym_n)
% [ii,k]=kgen(A,sym_n)
% this function generates the required connectivity information for the
% calculation of E_VSEPR.m

Natoms=size(A,1);
j1=uint16((1:Natoms/sym_n)');
[connec,connec2] = find(A(:,j1));
connec = uint16(connec);
j3 = next_node(j1,connec,3,A);
connec = reshape(connec,3,[])';                %% connectivity
connec2 = reshape(connec2,3,[])';                %% connectivity
j3=permute(reshape(j3,2,3,Natoms/sym_n),[3 2 1]);
k=zeros(Natoms/sym_n,3,3);
k(:,:,1)=connec;k(:,:,[2 3])=j3;

[ii,jj]=find(A);
ii=[(1:size(ii,1))' ii jj];
k1=kCompare(ii,[connec2(:,1) connec(:,1)]);
k2=kCompare(ii,[connec2(:,2) connec(:,2)]);
k3=kCompare(ii,[connec2(:,3) connec(:,3)]);
k11=kCompare(ii,[k(:,1,1) k(:,1,2)]);
k12=kCompare(ii,[k(:,1,1) k(:,1,3)]);
k21=kCompare(ii,[k(:,2,1) k(:,2,2)]);
k22=kCompare(ii,[k(:,2,1) k(:,2,3)]);
k31=kCompare(ii,[k(:,3,1) k(:,3,2)]);
k32=kCompare(ii,[k(:,3,1) k(:,3,3)]);
k=uint16([k1;k2;k3;k11;k12;k21;k22;k31;k32]);
end
function kk=kCompare(ii,k)
kk=zeros(length(k),1);
for j=1:length(k)
    kk(j)=ii(ii(:,2)==k(j,1)&ii(:,3)==k(j,2),1);
end

end
