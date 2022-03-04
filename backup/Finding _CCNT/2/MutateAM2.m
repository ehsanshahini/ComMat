function y=MutateAM2(x,mu,VarMax,VarMin,nVar)

    ni=VarMin;
    nj=VarMax; %hade balaye taghirat
    
    for i=1:nVar
       
        r1=rand;
        
        if r1<mu 
                while true
                    y(i)=randi([ni nj]);
                    if y(i)==x(i)
                        continue
                    end
                    break
                end
                if i==5
                    if x(5)==1
                       y(5)=2;
                    else
                       y(5)=1;                        
                    end
                end       
                if i==6
                    while true
                        y(i)=randi([1 9]);
                        if y(i)==x(i)
                            continue
                        end
                        break
                    end
                end
        else
            y(i)=x(i);
        end
    end
end