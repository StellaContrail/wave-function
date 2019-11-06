set grid
set xrange [-10:10]
set yrange [0:1]
set title "Probability density of BEC wave functions"
plot "data.txt" with lines title "final"
replot "data_initial.txt" with lines title "initial"
replot 0.5*x*x title "potential"