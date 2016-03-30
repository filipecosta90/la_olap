sparse = csvread('1_5_suppliers.csv');
row = sparse ( :, [1]); 
column = sparse ( :, [2]); 

n_rows = max(row);
n_columns = max(column);
c = zeros( n_rows , n_columns);

for i = 1 : n_rows
        c(row(i),column(i)) =  1;
end

non_zeros = 0;
zeros = 0;

for i = 1 : n_rows
    for j = 1 : n_columns
        if c(i,j) > 0 
                non_zeros = non_zeros + 1;
        else
                zeros = zeros +1 ;
        end
    end
end

percentage = non_zeros / zeros * 100;
spy(c)
