function alpha=LinearSearch(fun,alpha,s,X,E,g)
%http://www-fp.mcs.anl.gov/otc/Guide/OptWeb/continuous/unconstrained/linese
%arch.html

[Enew,gnew]=feval(fun,X+alpha*s);
gnew=gnew(:);
n=0;
para=0.9;
while Enew>E+0.001*alpha*g'*s || abs(gnew'*s)>para*abs(g'*s)
    alpha=alpha*0.8;
%     [Enew>E+0.001*alpha*g'*s abs(gnew'*s)>para*abs(g'*s)]
    Xnew=X+alpha*s;
    [Enew,gnew]=feval(fun,Xnew);
    gnew=gnew(:);
    n=n+1;
    if n>20
%         disp 'linesearch loop breaked'
        break
    end
end

end
