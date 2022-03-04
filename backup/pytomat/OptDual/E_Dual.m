function [E,E_gradient]=E_Dual(R,ii,n,sym_fun)

e0 = 1;
b0=1.43*sqrt(3);

R=feval(sym_fun,R,n);

Natoms=size(R,1);

uall=R(ii(:,2),:)-R(ii(:,1),:);
Ball=sqrt(sum(uall.^2,2));
stretching=0.5*e0*sum((Ball-b0).^2);
E = stretching;

if nargout > 1
    uall=uall./Ball;
    %% bond stretching term
%     E0_gradient=uall.*(Ball-b0);
    E0_gradient=NaN(Natoms,3);
    for k=1:Natoms
        ik=ii(:,2)==k;
        E0_gradient(k,:)=uall(ik,:)'*(Ball(ik)-b0);
    end
    E0_gradient= 2*e0*E0_gradient;
% keyboard
    E_gradient = E0_gradient;
    
    E_gradient = E_gradient(:);
end
end
