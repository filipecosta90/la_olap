SELECT l_returnflag, l_linestatus,  sum(l_quantity)
FROM LINEITEM_SAMPLE
GROUP BY l_returnflag, l_linestatus;
