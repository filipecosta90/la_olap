clear;

sparset = readtable('sparsity_9.csv','ReadVariableNames',true);
sparset.row = sparset.row + 1;
sparset.column = sparset.column +1;
n_columns = max(sparset.column);
values = 1:n_columns;

Q = sparse(sparset.row, sparset.column, values);



[m,n] = size(Q);

non_zeros = nnz (Q);

percentage = ( non_zeros / ( m * n ) ) * 100 ;


hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 0 960 960])
spy(Q);
set(gca,'YTickLabel',num2str(get(gca,'YTick')'));
set(gca,'XTickLabel',num2str(get(gca,'XTick')'));

str = sprintf('Sparsity pattern for Sparse Matrix based on GQuarks\nFor TPC-H 1 lineitem table, row number 9\nMatrix dimensions: ( %d x %d ), non-Zero Percentage: %.5f %%', m, n, percentage);
t = title(str);
set(t,'FontSize',14);
print('sparsity_9','-depsc');





