function A=conn(X)
eps=1e-5;
Natoms=size(X,1)/3;
[~,~,conn]=unique(round(X*1e5),'rows');
A=squeeze(sum(sum(reshape(abs(conn*(1./conn)'-1)<eps,[3 Natoms 3 Natoms]),1),3))==2;
end