
fid = fopen('lineitem.tbl');
fileID = fopen('15_base64_lineitem.csv','w');

pos = 15;
column = 0;
tline = fgets(fid);
while ischar(tline)
    C = strsplit(tline,'|');
   encoded = base64encode(char(C(pos)));
   finalpos = base64decode(encoded);
   row = sprintf('%d', finalpos);
   column = column +1 
   final_line  = sprintf('%s, %d, 1 ', row, column);
   fprintf(fileID, '%s\n', final_line);
   tline = fgets(fid);
end


fclose(fid);
fclose(fileID);

