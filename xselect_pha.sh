#!/bin/zsh
. ~/.zshrc
heainit

name_inevt=$1

xselect << EOF
MAXI_J1820p020
set datadir .
read event ${name_inevt}
y
extract spectrum
save spectrum
${name_inevt/.evt/}
exit
no
EOF
