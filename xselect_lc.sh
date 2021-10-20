#!/bin/zsh
. ~/.zshrc
heainit

name_inevt=$1
name_session=$2
ch_min=$3
ch_max=$4
ch_min_str=$5
ch_max_str=$6
binsize=$7

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
${name_inevt/.evt/_${ch_min_str}_${ch_max_str}}
exit
no
EOF
