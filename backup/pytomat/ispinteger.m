function TF=ispinteger(test,x)
% ispinteger(A) returns TRUE if (round(A)==A && isscalar(A) && A>0)
% returns TRUE.

test=test(:);

if nargin == 1
    if all(round(test)==test) && all(test>0)
        TF=logical(1);
    else
        TF=logical(0);
    end
else
    if all(round(test)==test) && all(test>0 | test==x) 
        TF=logical(1);
    else
        TF=logical(0);
    end
end 