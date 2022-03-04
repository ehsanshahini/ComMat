function R2=g_stripe(g,z)
if ~g
    R2=[];
    return
end
R=[0 0;1 0;0.5 1;0 0;0.5 1;-0.5 1];
R(:,1)=R(:,1)+0.5;
% R(R(:,2)==1,1)=R(R(:,2)==1,1)+0.5;
R=trans_expansion(0:z-1,R,1);
R2=[R;R(:,1) -R(:,2)];R2(:,2)=R2(:,2)+1;
R2=trans_expansion(2*(0:g/2-1),R2,2);
if mod(g,2)
    R2=[R2;R(:,1) -R(:,2)+g];
end
end
