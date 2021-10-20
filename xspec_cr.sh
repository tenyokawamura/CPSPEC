#!/bin/zsh
. ~/.zshrc
heainit

name_inpha=$1
name_outqdp=$2
low=$3
high=$4

xspec << EOF
data ${name_inpha}
iplot data
wd ${name_outqdp}
quit
exit
EOF
