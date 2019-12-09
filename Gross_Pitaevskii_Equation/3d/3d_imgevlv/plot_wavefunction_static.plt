set terminal png size 2400,800
set output "WF_static.png"
N    = 50
xmax = 10
set pm3d
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set xlabel "x (unit:Xs)"
set ylabel "y (unit:Xs)"
set zlabel rotate by 90

set multiplot layout 1,3
    set zrange  [0:1]
    set cbrange [0:1]
    set zlabel "Probability"
    set title "Initial trial wave function projected onto xy plane"
    splot "data_initial.txt" using 1:2:3 with pm3d title ""

    set title "Result wave function projected onto xy plane"
    splot "data_projection.txt" using 1:2:3 with pm3d title ""

    unset cbrange
    set zrange [-xmax:xmax]
    set zlabel "Z (Unit:Xs)"
    set title "External Potential cutted out"
    splot "data_potential_cutout.txt" every :::N*2::(N*3-1) using 1:2:3:($4*$4) with pm3d title ""
unset multiplot

unset output
set terminal wxt