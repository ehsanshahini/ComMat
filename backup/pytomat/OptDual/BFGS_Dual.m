function X = BFGS_Dual(R,A,n,sym_fun,flag)
% X = BFGS_Dual(R,A,n,sym_fun,flag)
% this function performs the optimization routine of a given graphitic
% molecule or a extended structure, with the VSEPR potential.
% Input R is the xyz coordinate matrix of the molecule. If the target
% molecular structure is periodic, the information of periodicity (e.g. the
% translation vector) should be included as R={R,l}, where l is the the
% periodicity information whose form depends on the symmetry function used.
% For aperiodic finite system R should simply be a matrix.
% A is the connectivity matrix of the molecule
% n is the symmetry number of the molecule. Depending on the symmetry
% function, n can be a scalar or a vector.
% sym_fun is the function handle of the symmetry function 
% flag='off' to not to display the opt. process, otherwise the process will
% be shown. flag is set to 'off' by default.

if nargin<5 % display opt. process
    flag='off';
end
sym_n=n(1);

[ii(:,1),ii(:,2)]=find(A);

X=R(1:end/sym_n,:);
[m1,n1]=size(X);
X=X(:);

%http://en.wikipedia.org/wiki/BFGS_method
invB=diag(ones(m1*n1,1)); % initial guess of Hessian matrix
[~,g]=E_Dual(reshape(X(1:m1*n1),[m1 n1]),ii,n,sym_fun);
alpha=0.1;
count=1;
if ~strcmp(flag,'off')
    disp(' ')
    disp('Broyden-Fletcher-Goldfarb-Shanno (BFGS) method with line search optimization started')
    disp(' ')
    disp('  nIters        f(x)           alpha            error')
end
E=E_Dual(reshape(X(1:m1*n1),[m1 n1]),ii,n,sym_fun);
if ~strcmp(flag,'off')
    disp(['       0    ' num2str(E,'%10.5e')])
end
while alpha>1e-5%||count<100

    s=-invB*g;
    alpha=LinearSearch(@(x) E_Dual(reshape(x(1:m1*n1),[m1 n1]),ii,n,sym_fun)...
    ,alpha*1.05,s,X,E,g);
    
    X=X+alpha*s;
    [E,newg]=E_Dual(reshape(X(1:m1*n1),[m1 n1]),ii,n,sym_fun);
    
    y=newg(:)-g(:);
    g=newg(:);
    error=norm(g);

    sy=s'*y;
    invBy=invB*y;
    sinvBy=s*invBy';
    invB=invB+(s*s')*(s'*y+y'*invBy)/sy^2-(sinvBy+sinvBy')/sy;
    
    if ~strcmp(flag,'off')
        format short;
        if ~mod(count,10),disp([repmat(' ',1,7-floor(log10(count))) num2str(count) '    ' num2str(E,'%10.5e') '    ' num2str(alpha,'%10.5e') '    ' num2str(error,'%10.5e')]),end
        if error>1e-5&&~mod(count,100),disp(' '),disp('  nIters        f(x)           alpha            error'),end;count=count+1;
    end
end

X=reshape(X(1:m1*n1),[m1 n1]);
X=feval(sym_fun,X,n);


end