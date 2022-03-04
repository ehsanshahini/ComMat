function [X,theta,u,t]= BFGS_helix(R,A,n,theta,u,flag)

if nargin<6
    flag='';
end

N=size(R,1);
% [m1,m2]=find(A);
% D=mean(sqrt(sum((R(m1,:)-R(m2,:)).^2,2)));
% R=R/D*1.43; % averaging the bond lengths
% % u=u/D*1.43;
% clear m1 m2 D
X=R(N/n+1:N/n*2,:);


j1=uint16((N/n+1:2*N/n)');
[connec,connec2] = find(A(:,j1));
connec = uint16(connec);
j3 = next_node(j1,connec,3,A);
connec = reshape(connec,3,[])';
connec2 = reshape(connec2,3,[])'+N/n;
j3=permute(reshape(j3,2,3,N/n),[3 2 1]);
k=zeros(N/n,3,3);
k(:,:,1)=connec;k(:,:,[2 3])=j3; 
% the first N-by-3 dimension of k is the nearest neighbors
% the second and third are the second nearest neighbors

[ii,jj]=find(A);
ii=[(1:3*N)' ii jj];
k1=kCompare(ii,[connec2(:,1) connec(:,1)]);
k2=kCompare(ii,[connec2(:,2) connec(:,2)]);
k3=kCompare(ii,[connec2(:,3) connec(:,3)]);
k11=kCompare(ii,[k(:,1,1) k(:,1,2)]);
k12=kCompare(ii,[k(:,1,1) k(:,1,3)]);
k21=kCompare(ii,[k(:,2,1) k(:,2,2)]);
k22=kCompare(ii,[k(:,2,1) k(:,2,3)]);
k31=kCompare(ii,[k(:,3,1) k(:,3,2)]);
k32=kCompare(ii,[k(:,3,1) k(:,3,3)]);
k=uint16([k1;k2;k3;k11;k12;k21;k22;k31;k32]);

tic
[m1,m2]=size(X);
X=[X(:);theta;u];

%http://en.wikipedia.org/wiki/BFGS_method
invB=diag(ones(m1*m2+2,1)); % initial guess of Hessian matrix
[E,g,sg]=helix_Bead_E(X,ii,k,n);

g=[g(:);sg(:)];
% error=1;
alpha=0.1;
count=1;
if ~strcmp(flag,'off')
    disp(' ')
    disp('Broyden-Fletcher-Goldfarb-Shanno (BFGS) method with line search optimization started')
    disp(' ')
    disp('  nIters        f(x)           alpha            error')
end
E=helix_Bead_E(X,ii,k,n);
if ~strcmp(flag,'off')
    disp(['       0    ' num2str(E,'%10.5e')])
end
while alpha>1e-5||count<100

    s=-invB*g;
    alpha=LinearSearch(@(x) helix_Bead_E(x,ii,k,n),alpha*1.05,s,X,E,g);
%     alpha=fminsearch(@(x) helix_Bead_E(X+x*s,ii,k,n),alpha);

    X=X+alpha*s;
    [E,ng,nsg]=helix_Bead_E(X,ii,k,n);
    
    y=[ng(:);nsg(:)]-g;
    g=[ng(:);nsg(:)];
    error=norm(g);
    
    sy=s'*y;
    invBy=invB*y;
    sinvBy=s*invBy';
    invB=invB+(s*s')*(sy+y'*invBy)/sy^2-(sinvBy+sinvBy')/sy;

    if ~strcmp(flag,'off')
        format short;
        if ~mod(count,10),disp([repmat(' ',1,7-floor(log10(count))) num2str(count) '    ' num2str(E,'%10.5e') '    ' num2str(alpha,'%10.5e') '    ' num2str(error,'%10.5e')]),end
        if error>1e-5&&~mod(count,100),disp(' '),disp('  nIters        f(x)           alpha            error'),end;count=count+1;
    end
end
t=toc;

theta=X(end-1);
u=X(end);
X=reshape(X(1:end-2),[m1 m2]);
end
function alpha=LinearSearch(fun,alpha,s,X,E,g)
%http://www-fp.mcs.anl.gov/otc/Guide/OptWeb/continuous/unconstrained/linese
%arch.html

[Enew,gnew,sgnew]=feval(fun,X+alpha*s);
gnew=[gnew(:);sgnew(:)];
n=0;
para=0.9;
while Enew>E+0.001*alpha*g'*s || abs(gnew'*s)>para*abs(g'*s)
    alpha=alpha*0.8;
%     [Enew>E+0.001*alpha*g'*s abs(gnew'*s)>para*abs(g'*s)]
    Xnew=X+alpha*s;
    [Enew,gnew,sgnew]=feval(fun,Xnew);
    gnew=[gnew(:);sgnew(:)];
    n=n+1;
    if n>20
%         disp 'linesearch loop breaked'
        break
    end
end

end
function kk=kCompare(ii,k)
kk=zeros(length(k),1);
for j=1:length(k)
    kk(j)=ii(ii(:,2)==k(j,2)&ii(:,3)==k(j,1),1);
end

end
function next = next_node(s,f,n,A)
s=repmat(s,1,n)';
s=s(:);

[t,null] = find(A(:,f));
t = reshape(t,3,[]);
s = repmat(s',3,1);
next=t(t~=s);

end
function r=range(t)
r=max(t)-min(t);
end