#!/bin/bash
#$1=city,$2=R,$3=lambda_limit,$4=seed
for i in 0 $3;do
./plot_cost_vs_lambda.sh log_lambda_vs_cost_${1}_R$2_S$4_LA  Cost_vs_Lambda_$1-R=$2-S=$4  R=$2-S=$4 $1-R=$2-S=$4 
done
