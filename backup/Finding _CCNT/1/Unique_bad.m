clear;
clc;

fileID = fopen('bad.txt','r');
numLines = 22;

A = [];
for k = 1:numLines
   mydata = fgetl(fileID);
   A=[A;mydata];
end
A=unique(A,'rows');
fclose all;

