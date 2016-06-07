BEGIN { FS=","; }
NR>=8 { print $1","$2","$13","$14 }
END { }
