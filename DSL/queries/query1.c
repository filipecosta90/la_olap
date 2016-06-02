BEGIN
bitmap returnFlag, lineStatus, shipdate_gt, shipdate_lt;
matrix quantity;

returnFlag = tbl_read('lineitem.tbl', 8);
lineStatus = tbl_read('lineitem.tbl', 9);
quantity = tbl_read('lineitem.tbl', 4);
shipdate_gt = tbl_filter('lineitem.tbl', 10, >=, '1998-08-28' );
shipdate_lt = tbl_filter('lineitem.tbl', 10, <=, '1998-12-01' );


bitmap projection, selection;
matrix aggregiation, intermediate_result, final_result;

vector bng;
bng = bang(6);

selection = shipdate_gt * shipdate_lt;
projection = returnFlag krao lineStatus;
aggregation = quantity * bng;
intermediate_result = projection * selection;
final_result = itermediate_result * aggregation;

tbl_write(final_result, 'final_result.tbl');
END


