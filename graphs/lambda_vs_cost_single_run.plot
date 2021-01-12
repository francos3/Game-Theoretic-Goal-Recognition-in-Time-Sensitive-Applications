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
set title font ",20"
#set title '{/Symbol l} Vs Cost,Shanghai,|D|=12'
set title '{/Symbol l} Vs Cost,100x100,|D|=8'

# Set label of x-axis
set xlabel '{/Symbol l}' font ",18"

# Set label of y-axis
set ylabel 'Cost' font ",18"

# Use a line graph
set style data line

#plot filename using 1:2:xticlabel(1) title columnheader with linespoints
plot filename using 1:2:xticlabel(1) notitle with linespoints

