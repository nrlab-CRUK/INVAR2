#!/bin/bash

DIR=$(dirname $0)

python3 $DIR/reorganiseInvarScores.py $1 .invarscoretemp.csv

Rscript $DIR/reorganiseInvarScores.R .invarscoretemp.csv $2 $3 $4

rm -f .invarscoretemp.csv

