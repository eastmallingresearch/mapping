#!/bin/sh


echo there are $# arguments
cd  arab_blast
echo $BLASTDIR/formatdb -i "$1" -n db -t "db" -p F
formatdb -i "$1" -n db -t "db" -p F


