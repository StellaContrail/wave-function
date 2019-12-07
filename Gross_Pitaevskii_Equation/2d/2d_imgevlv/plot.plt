set size ratio 1
set terminal jpeg size 1600,1600
set output "data.jpg"
set pm3d
set xrange  [-15:15]
set yrange  [-15:15]
set xlabel "X (unit:Xs)"
set ylabel "Y (unit:Xs)"
set zlabel rotate by 90
set multiplot layout 2,2
    set zrange  [0:1]
    set cbrange [0:0.01]
    set zlabel "Probability"
    set title "Result wave function"
    splot "data.txt" using 1:2:3 with pm3d title ""

    unset zrange
    unset cbrange
    set zlabel "Intensity"
    set title "External Potential"
    splot "data_potential.txt" using 1:2:3 with pm3d title ""

    set zrange [-1:1]
    set cbrange [-0.5:0.5]
    set style fill transparent solid 0.8
    set zlabel "Probability"
    set title "Real part of phased wave function"
    splot "data_phased.txt" using 1:2:4 with pm3d title ""
    set title "Imaginary part of phased wave function"
    splot "data_phased.txt" using 1:2:5 with pm3d title ""
    unset style
    unset cbrange
    unset zrange
unset multiplot
unset pm3d
unset output
unset size
set terminal wxt