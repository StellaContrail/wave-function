N    = 50
xmax = 10
set terminal wxt size 800,800
set size ratio 1
set xlabel "X (Unit:Xs)"
set ylabel "Y (Unit:Xs)"
set zlabel "Z (Unit:Xs)"
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set zrange [-xmax:xmax]
set zlabel rotate by 90
set multiplot layout 2,2
    set title "Potential at X=0"
    splot "data_potential.txt" every :::0::(N-1) using 1:2:3:($4*$4) with pm3d title ""
    
    set title "Potential at Y=0"
    splot "data_potential.txt" every :::N::(2*N-1) using 1:2:3:($4*$4) with pm3d title ""
    
    set title "Potential at Z=0"
    splot "data_potential.txt" every :::2*N::(3*N-1) using 1:2:3:($4*$4) with pm3d title ""
unset multiplot
unset zlabel
unset ylabel
unset xlabel