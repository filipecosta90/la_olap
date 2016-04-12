#!bin/bash

for file in *.tbl
do
  echo $file 
  sed -i -e "s/|/,/g" $file
  mv *.tbl-e $file.csv
done


