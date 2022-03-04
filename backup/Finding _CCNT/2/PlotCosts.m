function PlotCosts(pop1,pop2)

    Costs1=[pop1.Cost];
    Costs2=[pop2.Cost];
    close all;    
    hold on;
    plot3(Costs1(1,:),Costs1(2,:),Costs1(3,:),'r*','MarkerSize',8);
    plot3(Costs2(1,:),Costs2(2,:),Costs2(3,:),'bo','MarkerSize',8);
    xlabel('Radius');
    ylabel('Pitch Angle');
    zlabel('Diameter');
    grid on;
end