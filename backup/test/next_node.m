function next = next_node(s,f,n,A)
s=repmat(s,1,n)';
s=s(:);

[t,null] = find(A(:,f));
t = reshape(t,3,[]);
s = repmat(s',3,1);
next=t(t~=s);

end
