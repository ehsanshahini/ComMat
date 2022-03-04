clear;
clc;
close all;

load('pop1775final_sig_eps_E_SimNo.mat');

data = [];
for i=1:length(pop)
    Position = pop(i).Position;
    index = Position(1:4);
    type = Position(5);
    shift = Position(6);
    Parameters = helix(type,index,shift);
    data = [data;Parameters,-pop(i).Cost(1),-pop(i).Cost(2)]; % Parameter-Yield stress-Yield starin
end
