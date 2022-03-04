function R=sym_Cn(R,n)
% R=sym_Cn(R,n)
% this function performs the Cn rotational symmetry operation
% with z-axis the Cn rotation axis
% output R is n longer than input R

% warning('off','MATLAB:declareGlobalBeforeUse')
% global n
% warning('on','MATLAB:declareGlobalBeforeUse')

R=repmat(R,n,1);
for k=2:n
    Rot=RotatM([pi*2/n*(k-1) 0 0 1]);
    if size(R,2)==2,Rot=Rot(1:2,1:2);end
    R(end/n*(k-1)+1:end/n*k,:)=R(end/n*(k-1)+1:end/n*k,:)*Rot;
end
end