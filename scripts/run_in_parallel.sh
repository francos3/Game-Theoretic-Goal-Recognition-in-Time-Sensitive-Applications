#!/bin/bash
N=$1
R=$2
S=$3
LA1=$4
LA2=$5

for i in $(seq $4 1 $5)
do
  ./yen IL=0 N=$1 R=$2 S=$3 LA=$i optimize_budget min_cal MDD1=10000 MDD2=1 >& log_lambda_vs_cost_N${1}_R${2}_S${3}_LA_$i &
done
