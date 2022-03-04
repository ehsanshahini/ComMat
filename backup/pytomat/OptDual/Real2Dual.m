function [Rd,Ad,Nd,Tdual]=Real2Dual(R,T)
% [Rd,Ad,Nd,Tdual]=Real2Dual(R,T)
% this function returns the coordinates Rd, adjacency Ad, valency Nd, and
% the ring index matrix Td of a given molecule (R,T) in dual space.


NT=numel(T);
Nd=zeros(7,1);
for k=1:NT
   Nd(k)=size(T{k},1);
end
Nd=[0;cumsum(Nd)];

Rd=zeros(Nd(end),3);
for k=1:NT
    Rd(Nd(k)+1:Nd(k+1),:)=squeeze(mean(reshape(R(T{k}',:),k+2,[],3),1));
end
N=Nd(end);

Ad=false(N,N);
Tdual=NaN(size(R,1),3);
for k1=1:size(R,1)
    m=[];
    for k2=1:NT
        if ~isempty(T{k2})
            m=[m;find(any(T{k2}==k1,2))+Nd(k2)];
        end
    end
    
    switch numel(m)
        case 2
            Ad(m(1),m(2))=true;
        case 3
            Ad(m(1),m(2))=true;Ad(m(2),m(3))=true;Ad(m(3),m(1))=true;
            Tdual(k1,:)=m;
%         warning('AAA:BBB',['atom No.' num2str(k1) ' found not to be shared by three faces'])
    end
    
end
% keyboard
Ad=Ad|Ad';

% m=1;
% Tdual=cell(7,1);
% return
% for kk=1:7
%     if Nd(kk+1)==Nd(kk)
%         continue
%     end
%     Ring_xyz=zeros(Nd(kk+1)-Nd(kk),kk+2,3);
%     for k=1:(Nd(kk+1)-Nd(kk))
%         Tdual{kk}(k,:)=find(Ad(m,:));
%         Ring_xyz(k,:,:)=bsxfun(@minus,Rd(Tdual{kk}(k,:),:),Rd(m,:));
%         Ring_xyz(k,:,:)=bsxfun(@rdivide,Ring_xyz(k,:,:),dnorm(Ring_xyz(k,:,:),3));
%         for l=1:kk+1
%             [null,ii]=max(sum(bsxfun(@times,Ring_xyz(k,l,:),Ring_xyz(k,l+1:end,:)),3));
%             Tdual{kk}(k,[l+1 ii+l])=Tdual{kk}(k,[ii+l l+1]);
%             Ring_xyz(k,[l+1 ii+l],:)=Ring_xyz(k,[ii+l l+1],:);
%         end
%         m=m+1;
%     end
% end

end