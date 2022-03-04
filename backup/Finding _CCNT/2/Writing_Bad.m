function Writing_Bad(x)

fileID1=fopen('bad.txt','r');
A=fscanf(fileID1,'%f');

fileID2=fopen('bad.txt','a');

k=sprintf('%1d',x);
iwant=str2double(k);

if ~ismember(iwant,A)
    fprintf(fileID2,'%d',x);
    fprintf(fileID2,'\n');
end

fclose(fileID1);
fclose(fileID2);


end
