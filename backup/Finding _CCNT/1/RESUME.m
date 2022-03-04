

for i=1:length(poptotal)        %end of poptotal
    pop(i).Position=poptotal(i).Position;
    pop(i).Cost(1,1)=-poptotal(i).YieldStress;  
    pop(i).Cost(2,1)=-poptotal(i).Strainend;
end



% Non-Dominated Sorting
    [pop F]=NonDominatedSorting(pop);

    % Calculate Crowding Distance
    pop=CalcCrowdingDistance(pop,F);

    % Sort Population
    [pop F]=SortPopulation(pop);
    
        % Store F1
    F1=pop(F{1});
    
%         figure(1);
%     PlotCosts(F1);
    
    
    
    
    
    