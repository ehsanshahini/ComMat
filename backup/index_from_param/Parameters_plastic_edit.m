clear;
clc;
close all;

load('poptotal_sigY_T_5617_plastic.mat');
load('data_Plastic.mat')

data_Plastic = [];
for i=1:length(poptotal)
    Position = poptotal(i).Position;
    index = Position(1:4);
    type = Position(5);
    shift = Position(6);
    Parameters = helix(type,index,shift);
    data_Plastic = [data_Plastic;Parameters,poptotal(i).YieldStress,poptotal(i).Strainend,poptotal(i).Toughness];
end