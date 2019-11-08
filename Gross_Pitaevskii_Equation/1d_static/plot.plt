set grid
set xrange [-5:5]
set yrange [0:1]
sigma = 0.05
set title "Probability density of BEC wave functions"
plot "data_initial.txt" with lines title "initial"
replot "data.txt" with lines title "final"
replot 0.5*x*x+exp(-0.5*x*x/(sigma*sigma)) title "potential"