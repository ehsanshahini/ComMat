function T=ring_ind_c(A,M)

% [c1,c2]=find(A);
% N=size(A,1);


if nargin<2
    M=7;
end
A=sparse(A);
T=ring_ind_for(A,M);
T=non_elementary(A,T,M);
if M>6 % check for 6-member rings running across any (3,0) nanotube
    T=check30CNT(T);
end
if M>7 % check for 8-member rings running across any (4,0) nanotube
    T=check40CNT(T);
end

end


function T=non_elementary(A,T,M)

An=double(A);
D=uint16(full(A));
for i=2:M
    An=An*A;
    D(An&~D)=i;
end

for i=1:M-2
    non_elementary=[];
    for j=1:size(T{i},1)
        t=0;
        for k=1:i+2
            for l=k+1:i+2
                if D(T{i}(j,k),T{i}(j,l))<min(l-k,i+2-l+k)
                    non_elementary=[non_elementary j];
                    t=1;
                    break
                end
            end
            if t==1;break,end;
        end
    end
    T{i}(non_elementary,:)=[];
end


end

function T=check40CNT(T)
T6=T{4};
T8=T{6};

n8=size(T8,1);
ind=false(n8,1);
for k8=1:n8
    ind(k8)=any(sum(sum(bsxfun(@minus,permute(T8(k8,:),[1 3 2]),T6)==0,3),2)>2);
end
T{6}(ind,:)=[];
end
function T=check30CNT(T)
T6=T{4};

n6=size(T6,1);
ind=false(n6,1);
for k6=1:n6
    ind(k6)=sum(sum(sum(bsxfun(@minus,permute(T6(k6,:),[1 3 2]),T6)==0,3),2)>2)>3;
end
T{4}(ind,:)=[];
end