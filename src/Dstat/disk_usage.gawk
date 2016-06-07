BEGIN { FS=","; }
NR>=8 { print $1","$2","$11","$12 }
END { }
