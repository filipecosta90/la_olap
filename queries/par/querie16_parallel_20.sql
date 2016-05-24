set max_parallel_degree = 20;
\timing on
SELECT l_returnflag, l_linestatus,  sum(l_quantity)
FROM LINEITEM_16
WHERE l_shipdate >= '1998-08-28' AND l_shipdate <= '1998-12-01' 
GROUP BY l_returnflag, l_linestatus;
\timing off
