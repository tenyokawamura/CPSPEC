#!/bin/zsh
. ~/.zshrc
heainit

name_inevt=$1
name_session=$2

xselect << EOF
${name_session}
set datadir .
read event ${name_inevt}
y
extract spectrum
save spectrum
${name_inevt/.evt/}
exit
no
EOF
