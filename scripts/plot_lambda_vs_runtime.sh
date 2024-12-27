#!/bin/bash
echo "lambda,$2" > results_$2.csv
grep overall $1* | cut -d ',' -f 35,14 | awk -F "," '{ print $2 "," $1}' | sed 's/ //g' | sort -n >> results_$2.csv
gnuplot -e "filename='results_$2.csv'" lambda_vs_runtime_single_run.plot >& $2.pdf
