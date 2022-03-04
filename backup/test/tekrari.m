clear;
clc;
close all;

pop1775 = load('pop1775_Elastic_small_cost=stress_strain.mat');
pop3031 = load('pop3031_Elastic_mix_cost=stress_strain.mat');

count = 0;
k = 1;
for i=1:length(pop1775.pop)
    count = 0;
    for j=1:length(pop3031.pop)
        
        if isequal(pop1775.pop(i).Position,pop3031.pop(j).Position)
            count = 1;
        end
    end
    if count == 0
        a(k)=pop1775.pop(i);
        k = k +1;
    end
    
end


% count = 0;
% for i=1:length(pop1775.pop)
%     for j=1:length(pop3031.pop)
%         
%         if isequal(pop1775.pop(i).Position,pop3031.pop(j).Position)
%             count = count + 1;
%         end
%     end
%     
% end
