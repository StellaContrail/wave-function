set grid
set xrange [-10:10]
set yrange [0:1]
set title "Probability density of BEC wave functions"

plot "data_initial.txt" using 1:2 with lines title "initial"
#replot "data_shifted.txt" using 1:4 with lines title "shifted"
replot "data.txt" using 1:2 with lines title "final"
replot "data_potential.txt" using 1:2 with lines title "potential"

set ylabel "Intensity"
set xlabel "x (unit:Xs)"