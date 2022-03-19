function v=trans_expansion(T,v,j)
temp=repmat(T,size(v,1),1);temp=temp(:);
v=repmat(v,length(T),1);v(:,j)=v(:,j)+temp;