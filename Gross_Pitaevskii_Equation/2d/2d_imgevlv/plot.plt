set pm3d
set multiplot layout 1,3
set xrange [-15:15]
set yrange [-15:15]
set zrange [0:1]
#set cbrange [0:0.3]
set xlabel "x (unit:Xs)"
set ylabel "y (unit:Xs)"
set zlabel "Intensity"
set zlabel rotate by 90
set title "Initial trial wave function"
splot "data_initial.txt" using 1:2:3 with pm3d title "initial"
set title "Result wave function"
splot "data.txt" using 1:2:3 with pm3d title "final"
#set title "Result wave function half of whose phase shifted by PI"
#splot "data_shifted.txt" using 1:2:4 with pm3d title "shifted (real)"
#splot "data_shifted.txt" using 1:2:5 with pm3d title "shifted (imag)"
unset zrange
#unset cbrange
set title "External Potential"
splot "data_pot.txt" using 1:2:3 with pm3d title "potential"
unset multiplot