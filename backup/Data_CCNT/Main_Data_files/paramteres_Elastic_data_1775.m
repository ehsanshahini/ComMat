clear;
clc;
close all;

load('a.mat');

data = [];
for i=1:length(a)
    Position = a(i).Position;
    index = Position(1:4);
    type = Position(5);
    shift = Position(6);
    Parameters = helix(type,index,shift);
    data = [data;Parameters,-a(i).Cost(1),-a(i).Cost(2)]; % Parameter-Yield stress-Yield starin
end
