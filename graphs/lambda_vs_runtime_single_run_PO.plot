# Output W3C Scalable Vector Graphics
#set terminal jpeg
#set terminal svg
#set terminal dumb
#set terminal latex
set terminal pdf
set grid

#legend placement
set key right top
set key vertical
set key font ",10"


# Read comma-delimited data from file
set datafile separator comma

# Set graph title
set title font ",18"
set title '{/Symbol l} Vs Runtime,100x100,|D|=8'
#set title '{/Symbol l} Vs Runtime,Shanghai,|D|=12'

# Set label of x-axis
set xlabel '{/Symbol l}' font ",14"

# Set label of y-axis
set ylabel 'RunTime(ms)' font ",14"

# Use a line graph
#set style data line
#set linestyle  1 linetype  1 linewidth 8
#set linestyle  2 linetype  2 linewidth 8
#set linestyle  3 linetype  3 linewidth 8
#set linestyle  4 linetype  4 linewidth 8
#set linestyle  5 linetype  5 linewidth 8
#set linestyle  6 linetype  6 linewidth 8
#set linestyle  7 linetype  7 linewidth 8
#set linestyle  8 linetype  8 linewidth 8
#set linestyle  9 linetype  9 linewidth 8
#set linestyle 10 linetype 10 linewidth 8

#plot filename using 1:2:xticlabel(1) title columnheader with linespoints
#plot filename using 1:2:xticlabel(1) notitle with linespoints
#plot for[i in "PO=45 PO=50 PO=60 PO=70 PO=80 PO=90" ] filename using 1:i:xticlabel(1) title columnheader with linespoints ps 0.6
#plot for[i in "PO=5 PO=10 PO=45 PO=50 PO=60 PO=70 PO=80 PO=90 PO=100" ] filename using 1:i:xticlabel(1) title columnheader with linespoints ps 0.4
plot for[i in "PO=45 PO=50 PO=60 PO=70 PO=80 PO=90 PO=100" ] filename using 1:i:xticlabel(1) title columnheader with linespoints ps 0.4
#plot for[i in "PO=70 PO=75 PO=80 PO=90 PO=95 PO=100" ] filename using 1:i:xticlabel(1) title columnheader with linespoints ps 0.4

