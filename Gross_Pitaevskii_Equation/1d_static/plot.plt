set grid
set xrange [-10:10]
set yrange [0:1]
sigma = 0.1
set title "Probability density of BEC wave functions"
plot "data_initial.txt" using 1:4 with lines title "initial"
replot "data_shifted.txt" using 1:4 with lines title "shifted"
replot "data.txt" using 1:4 with lines title "final"
replot 0.5*x*x+100*exp(-0.5*x*x/(sigma*sigma)) title "potential"