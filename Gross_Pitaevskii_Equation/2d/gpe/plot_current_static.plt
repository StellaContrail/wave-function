set terminal wxt size 800,800
load "gnuplot_vars.txt"
# SETTINGS
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set title  "FLUX"
set xlabel "X/a0"
set ylabel "Y/a0"
set palette rgbformulae 22,13,10
stats fn_wavefunction_imaginary_result using 3 nooutput
set zrange [STATS_min:STATS_max]
set cbrange [STATS_min:STATS_max]
stats fn_current_imag_result using 3:4 nooutput
SCALE = 1.0 / sqrt(STATS_max_x**2 + STATS_max_y**2)
set pm3d map
splot fn_wavefunction_imaginary_result using 1:2:3 title "", fn_current_imag_result using 1:2:(0):($3*SCALE):($4*SCALE):(0) title "" with vectors linecolor rgb "#000000"
