function y=MutateAM(x,mu,VarMax)

    nVar=numel(x);
    nj=VarMax; %hade balaye taghirat
    
    for i=1:nVar
       
        r1=rand;
        
        if r1<mu  %0.2 is same as mu in nsga2
            r2=rand;         
            if r2<0.5
                y(i)=max(1,x(i)-1);
                if i==5
                    if x(5)==1
                       y(5)=2;
                    else
                       y(5)=1;                        
                    end
                end       
            else 
            y(i)=min(nj,x(i)+1);
                if i==6
                    y(i)=min(9,x(i)+1);
                end
                if i==5
                    if x(5)==1
                       y(5)=2;
                    else
                       y(5)=1;                        
                    end
                end      
            end 
        else
            y(i)=x(i);
        end
    end
end