# SETTINGS
load "gnuplot_vars.txt"
set xrange [-xmax:xmax]
set yrange [0:1]
set xlabel "X/a0"
set ylabel "|{/Symbol y}|^2"
set xtics 1
set title "Density"
plot fn_wavefunction_imaginary_result using 1:2 title "WF" with lines
replot fn_potential_imaginary using 1:2 title "POTENTIAL" with lines