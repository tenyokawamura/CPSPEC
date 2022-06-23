#!/bin/zsh
. ~/.zshrc
heainit

name_inevt=$1
name_outpha=$2
name_session=$3

xselect << EOF
${name_session}
set datadir .
read event ${name_inevt}
y
extract spectrum
save spectrum
${name_outpha}
exit
no
EOF
