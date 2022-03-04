clear;
clc;
close all;

load('a.mat');

data_el_from_pl = [];
for i=1:length(a)
    Position = a(i).Position;
    index = Position(1:4);
    type = Position(5);
    shift = Position(6);
    Parameters = helix(type,index,shift);
    data_el_from_pl(i,:) = [Parameters,a(i).YieldStress,a(i).Strainend,a(i).Toughness,a(i).YieldStrain,a(i).YoungModulus];
end

% load('pop3031_sigma_eps.mat');
% data_Elastic = [];
% for i=2013:length(pop)
%     Position = pop(i).Position;
%     index = Position(1:4);
%     type = Position(5);
%     shift = Position(6);
%     Parameters = helix(type,index,shift);
%     data_Elastic = [data_Elastic;Parameters,-pop(i).Cost(1),-pop(i).Cost(2)]; % Parameter-Yield stress-Yield starin
% end