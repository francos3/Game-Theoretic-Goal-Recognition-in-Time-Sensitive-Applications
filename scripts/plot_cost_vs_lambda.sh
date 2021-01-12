#!/bin/bash
echo "$3,lambda," > results_$2.csv
grep overall $1* | cut -d ',' -f 35,10 | awk -F "," '{ print $2 "," $1}' | sed 's/ //g' | sort -n >> results_$2.csv
gnuplot -e "filename='results_$2.csv'" cost_vs_lambda_single_run.plot >& $2.pdf
