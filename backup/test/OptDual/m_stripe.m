function R=m_stripe(t1,t2,n)
% R=m_stripe(t1,t2,n)
% this function produces trapezoidal patch of graphene with larger base
% t2 and smaller base t1, thus t2>t1, n fot the rotational symmetry number, 6 is set
% for default which corresponds to the original honeycomb lattice

if nargin<3
    n=6;
end

err=1e-3;
R=[0 0;1 0;0.5 1;1 0;0.5 1;1.5 1];
R=[R;R(:,1) 2-R(:,2)];
R=trans_expansion(0:t2-1,R,1);
R=trans_expansion(0:2:2*(t2-t1-1),R,2);
R(:,2)=R(:,2)*sqrt(3)/2;

x=reshape(R(:,1),3,[]);
y=reshape(R(:,2),3,[]);
ind=reshape(repmat(all(y>sqrt(3)*x-err|y>(t2-t1)*sqrt(3)/2-err|y>-sqrt(3)*x+sqrt(3)*t2-err,1),3,1),[],1);
R(ind,:)=[];
R(:,1)=R(:,1)-t2/2;
R(:,2)=R(:,2)/sqrt(3)*cot(pi/n)-t2*cot(pi/n)/2;
end