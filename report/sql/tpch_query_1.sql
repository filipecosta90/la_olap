SELECT l_returnflag, l_linestatus,  sum(l_quantity)
FROM LINEITEM_SAMPLE
WHERE l_shipdate >= "1996-04-12" AND l_shipdate <= "1997-01-28"
GROUP BY l_returnflag, l_linestatus;
