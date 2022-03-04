clear;
clc;
close all;

popplastic = load('poptotal_sigY_T_5617_plastic.mat');
popelastic = load('pop_el_tot.mat');
load('data_Plastic.mat');


count = 0;
k = 1;
for i=1:length(popplastic.poptotal)
    count = 0;
    for j=1:length(popelastic.pop_el_tot)
        
        if isequal(popplastic.poptotal(i).Position,popelastic.pop_el_tot(j).Position)
            count = 1;
        end
    end
    if count == 0
        a(k,:)=data_Plastic(i,:);
        k = k +1;
    end
end
    
