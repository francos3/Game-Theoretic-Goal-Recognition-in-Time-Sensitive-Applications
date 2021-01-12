# Output W3C Scalable Vector Graphics
#set terminal jpeg
#set terminal svg
#set terminal dumb
#set terminal latex
set terminal pdf

#legend placement
set key top
set key horizontal
set key inside
set key font ",14"
set key maxrows 2

# Read comma-delimited data from file
set datafile separator comma

# Set graph title
set title font ",20"
set title 'Runtime as a function of grid size and |D|.'

# Set label of x-axis
set xlabel "Destinations" font ",18"

# Set label of y-axis
set ylabel "average time(secs)" font ",18"


# Use a line graph
set style data line

plot for [i=2:5] 'merged_results_runtime.csv' using 1:i title columnheader with linespoints
