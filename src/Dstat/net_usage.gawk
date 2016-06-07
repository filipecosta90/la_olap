BEGIN { FS=","; }
NR>=8 { print $1","$2","$(NF-8)","$(NF-7) }
END { }
