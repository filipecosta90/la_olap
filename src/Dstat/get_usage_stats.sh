#!bin/bash
CURRENT=`pwd`
BASENAME=`basename "$CURRENT"`

mkdir __memory_usage
mkdir __system_usage
mkdir __disk_usage
mkdir __net_usage
mkdir __cpu_usage

cd __csv/

for file in *.csv
do
  echo $file
  gawk -f ../memory_usage.gawk $file > ../"MEMORY_"$file
  gawk -f ../system_usage.gawk $file > ../"SYSTEM_"$file
  gawk -f ../disk_usage.gawk $file > ../"DISK_"$file
  gawk -f ../net_usage.gawk $file > ../"NET_"$file
  gawk -f ../cpu_usage.gawk $file > ../"CPU_"$file
done

cd ..

mv MEMORY_* __memory_usage
mv SYSTEM_* __system_usage
mv DISK_* __disk_usage
mv NET_* __net_usage
mv CPU_* __cpu_usage

