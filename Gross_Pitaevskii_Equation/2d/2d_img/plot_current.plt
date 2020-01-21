# Gnuplot script : Plot probability current
set terminal jpeg size 1600,1600
set size ratio 1
# SETTINGS
xmax        = 8     # BOUNDARY OF X
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set xlabel "X (Unit:Xs)"
set ylabel "Y (Unit:Xs)"
set view 0, 0, 1, 1
set palette rgbformulae 22,13,10
set pm3d map
set title sprintf("Probability current")
set output "current_initial.jpg"
splot "data_initial.txt" using 1:2:3 title "", "data_initial_flux.txt" using 1:2:(0):($3*5):($4*5):(0) title "" with vectors linecolor rgb "#000000"
set output "current_final.jpg"
splot "data.txt" using 1:2:3 title "", "data_flux.txt" using 1:2:(0):($3*5):($4*5):(0) title "" with vectors linecolor rgb "#000000"
unset pm3d
unset output
set palette rgbformulae 7,5,15
unset view
set terminal wxt