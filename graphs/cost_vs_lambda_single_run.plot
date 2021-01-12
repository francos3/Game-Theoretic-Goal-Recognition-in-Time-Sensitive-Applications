# Output W3C Scalable Vector Graphics
#set terminal jpeg
#set terminal svg
#set terminal dumb
#set terminal latex
set terminal pdf
set grid

#legend placement
set key right top

# Read comma-delimited data from file
set datafile separator comma

# Set graph title
#set title 'Cost Vs Lambda'
set title '{/Symbol l} Vs Cost,100x100,|D|=12'

# Set label of x-axis
set xlabel 'Cost'

# Set label of y-axis
set ylabel '{/Symbol l}'

# Use a line graph
set style data line

plot filename using 2:1 title columnheader with linespoints
plot filename using 2:1 notitle with linespoints

