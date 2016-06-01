clear;
sparse = csvread('__sparsity_csv/1_2_lineitem.csv');
row = sparse ( :, [1]); 
column = sparse ( :, [2]); 

n_rows = max(row);
n_columns = max(column);
c = zeros( n_rows , n_columns);

for i = 1 : length(row)
        c(row(i),column(i)) =  1;
end

non_zeros = 0;
zeros = 0;

for i = 1 : length(row)
    for j = 1 : length(column)
        if c(i,j) > 0 
                non_zeros = non_zeros + 1;
        else
                zeros = zeros +1 ;
        end
    end
end

percentage = non_zeros / zeros * 100;

hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 0 960 960])
spy(c,20);
title('Sparsity pattern for Sparse Matrix based on GNU QUARKS, TPC-H 1, suppliers table, index row 1, data row 6')

