set terminal wxt size 800,800
set pm3d map
load "gnuplot_vars.txt"
# SETTINGS
M           = N + 1
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set title "Density distribution of superfluid"
set xlabel "X/a0"
set ylabel "Y/a0"
set xtics 1
set ytics 1
stats fn_wavefunction_imaginary_result u 3 nooutput
set zrange  [0:STATS_max]
set cbrange [0:STATS_max]
splot fn_wavefunction_imaginary_result using 1:2:3:3 title ""
unset output
set terminal wxt
unset pm3d