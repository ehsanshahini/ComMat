clear;
clc;
close all;

popplastic = load('poptotal_sigY_T_5617_plastic.mat');
popelastic = load('pop_el_tot.mat');


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
        a(k)=popplastic.poptotal(i);
        k = k +1;
    end
    
end

% 
% count = 0;
% for i=1:length(pop5617.poptotal)
%     for j=1:length(pop3031.pop)
%         
%         if isequal(pop5617.poptotal(i).Position,pop3031.pop(j).Position)
%             count = count + 1;
%             display(pop5617.poptotal(i).Position);
%         end
%     end
%     
% end
