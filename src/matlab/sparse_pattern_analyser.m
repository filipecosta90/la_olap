load __sparsity_csv/1_2_lineitem.csv
H = spconvert(X1_2_lineitem);
load base64testing/15_base64_lineitem.csv
B = spconvert(X15_base64_lineitem);

hFig = figure(1);
set(gcf,'PaperPositionMode','auto')
set(hFig, 'Position', [0 0 960 960])
spy(H,'r');
hold on;
spy(B,'b');
hold off;

[m,n] = size(H);
[mB,nB] = size(B);

non_zeros = nnz (H);
non_zerosB = nnz (B);

percentage = ( non_zeros / ( m * n ) ) * 100 ;
percentageB = ( non_zerosB / ( mB * nB ) ) * 100 ;

set(gca,'YTickLabel',num2str(get(gca,'YTick')'));
set(gca,'XTickLabel',num2str(get(gca,'XTick')'));


str = sprintf('Sparsity pattern for Sparse Matrix based on GNU QUARKS vs Base64 Encoding,\nTPC-H 1 lineitem table, projection row 15, Non-Zero Percentage: %.5f %% GNU QUARKS vs %.15f %% Base64 Encoding,\nmatrices dimensions: ( %d x %d ) GNU QUARKS vs ( %d x %d ) Base64 Encoding',percentage, percentageB, m, n, mB, nB);
title(str);

