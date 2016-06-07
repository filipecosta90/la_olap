BEGIN { FS=","; }
NR>=8 { print $1","$2","$5","$8 }
END { }
