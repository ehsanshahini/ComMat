function R=m_stripe_tri(t1,t2)
% R=m_stripe_tri(t1,t2)
% version for the construction of triangular necks of m_stripe.m
R=[0 0;1 0;0.5 1;1 0;0.5 1;1.5 1];
R=[R;R(:,1) 2-R(:,2)];
R=trans_expansion(0:t2*3+t1,R,1);
R=trans_expansion(0:2:t2*3+t1,R,2);
R(:,2)=R(:,2)*sqrt(3)/2;
x=mean(reshape(R(:,1),3,[]),1);
y=mean(reshape(R(:,2),3,[]),1);
err=1e-5;
ind=reshape(repmat(y>x/sqrt(3)+err|y>t2*sqrt(3)/2|y>(-x+t1+t2*3)/sqrt(3)-err,3,1),[],1);
R(ind,:)=[];
R(:,1)=R(:,1)-(t2*3+t1)/2;
R(:,2)=R(:,2)-(t2*3+t1)*sqrt(3)/6;
end