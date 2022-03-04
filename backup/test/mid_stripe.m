function R=mid_stripe(z,m)
% R=mid_stripe(z,m)
% the horizontally squeezed version of m_stripe.m
% this function is useful for the construction of graphitic structures
% through parametric way instead of the direct real space method

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
        R(abs(R(:,2)-k)<eps&abs(R(:,1)-l)<eps,1)=R(abs(R(:,2)-k)<eps&abs(R(:,1)-l)<eps,1)-k*l/(z+k);
    end
end
R(:,1)=R(:,1)/2;
end
