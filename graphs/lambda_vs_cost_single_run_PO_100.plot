# Output W3C Scalable Vector Graphics
#set terminal jpeg
#set terminal svg
#set terminal dumb
#set terminal latex
set terminal pdf
set grid

#legend placement
set key top right
set key vertical
set key font ",16"
#set key outside
#set key maxrows 1

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
#set style data line
#set style line 1 lw 4 ps 2 pt 5 lc rgb "red"
#>    set style line 2 lw 4 ps 2 pt 7 lc rgb "blue"
#>    set style line 3 lw 4 ps 2 pt 9 lc rgb "orange"
#>    set style line 4 lw 4 ps 2 pt 11 lc rgb "web-green"
#>    set style line 5 lw 4 ps 2 pt 13 lc rgb "purple"
#>    set style line 6 lw 4 ps 2 pt 15 lc rgb "brown"

#plot filename using 1:2:xticlabel(1) title columnheader with linespoints
#plot filename using 1:2:xticlabel(1) notitle with linespoints
#plot for[i in "PO=45 PO=50 PO=60 PO=70 PO=80 PO=90" ] filename using 1:i:xticlabel(1) title columnheader with linespoints ps 0.6
#plot for[i in "PO=70 PO=75 PO=80 PO=85 PO=90 PO=95 PO=100" ] filename using 1:i:xticlabel(1) title columnheader with linespoints ps 0.6
#plot for[i in "PO=5 PO=10 PO=45 PO=50 PO=60 PO=70 PO=80 PO=90 PO=100" ] filename using 1:i:xticlabel(1) title columnheader with linespoints ps 0.5
plot for[i in "PO=50 PO=60 PO=70 PO=80 PO=90 PO=100" ] filename using 1:i:xticlabel(1) title columnheader with linespoints ps 0.4

