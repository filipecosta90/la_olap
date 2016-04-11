clear;
load __sparsity_csv/1_9_lineitem.csv
H = spconvert(X1_9_lineitem);

hFig = figure(1);
%set(gcf,'PaperPositionMode','auto')
%set(hFig, 'Position', [0 0 960 960])
spy(H);
[m,n] = size(H);
non_zeros = nnz (H);
percentage = ( non_zeros / ( m * n ) ) * 100 ;
set(gca,'YTickLabel',num2str(get(gca,'YTick')'));
set(gca,'XTickLabel',num2str(get(gca,'XTick')'));

str = sprintf('Sparsity pattern for Sparse Matrix based on GNU QUARKS, TPC-H 1\n lineitem table, index row 1, data row 9, Non-Zero Percentage: %.5f %%,\n dimensions matrix (%d,%d)',percentage, m, n);
title(str);

