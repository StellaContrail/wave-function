set grid
set xrange [-5:5]
set yrange [0:1]
sigma = 0.05
set sample 1000
set title "Probability density of BEC wave functions"
plot "data_initial.txt" using 1:4 with lines title "initial"
replot "data.txt" using 1:4 with lines title "final"
replot 0.5*x*x title "potential"