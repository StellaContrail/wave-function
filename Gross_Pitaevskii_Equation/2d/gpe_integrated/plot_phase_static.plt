set terminal wxt size 800,800
load "gnuplot_vars.txt"
# SETTINGS
set xtics 1
set ytics 1
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set title "PHASE DISTRIBUTION"
set xlabel "X/a0"
set ylabel "Y/a0"
stats fn_phase_distribution_imag_result using 4 nooutput
set zrange  [STATS_min:STATS_max]
set cbrange [STATS_min:STATS_max]
set pm3d map
splot fn_phase_distribution_imag_result using 1:2:4 title ""
unset pm3d