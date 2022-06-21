#!/bin/zsh
. ~/.zshrc
heainit

name_inevt=$1
name_outlc=$2
name_session=$3
ch_min=$4
ch_max=$5
binsize=$6

xselect << EOF
${name_session}
set datadir .
read event ${name_inevt}
y
set binsize ${binsize}
set phaname PI
filter pha_cutoff
${ch_min}
${ch_max}
extract curve
save curve
${name_outlc}
exit
no
EOF
