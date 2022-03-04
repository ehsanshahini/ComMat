function R=g0_tri(z)
% close all
% clear all
% z=3;

R=[0 0;1 0;0.5 1;1 0;0.5 1;1.5 1];
R=[R;R(:,1) 2-R(:,2)];
R=trans_expansion(0:z-1,R,1);
R=trans_expansion(0:2:z-1,R,2);
x=mean(reshape(R(:,1),3,[]));
y=mean(reshape(R(:,2),3,[]));
ind=reshape(repmat(y>2*x|y>-2*x+2*z,3,1),[],1);
R(ind,:)=[];
R=[R(:,2) R(:,1)/z-0.5];
% dplot(R)
end