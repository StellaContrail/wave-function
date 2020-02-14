# Gnuplot script : Plot result wave function density, real, imaginary part, and potential form
set size ratio 1
set terminal jpeg size 1600,1600
set output "data_initial.jpg"
set pm3d
xmax = 3
set xrange  [-xmax:xmax]
set yrange  [-xmax:xmax]
set xlabel "X (unit:Xs)"
set ylabel "Y (unit:Xs)"
set zlabel rotate by 90
set multiplot layout 2,2
    set zlabel "Probability"
    set title "Initial wave function"
    splot "data_initial.txt" using 1:2:3 with pm3d title ""

    set style fill transparent solid 0.8
    set zlabel "Intensity"
    set title "External Potential"
    splot "data_potential.txt" using 1:2:4 with pm3d title ""

    set palette defined ( -1 '#000030', 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#90ff70',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
    set zlabel "Probability"
    set title "Real part of initial wave function"
    splot "data_initial.txt" using 1:2:4 with pm3d title ""
    set title "Imaginary part of initial wave function"
    splot "data_initial.txt" using 1:2:5 with pm3d title ""
    unset style
    set palette rgbformulae 7,5,15
unset multiplot
unset pm3d
unset output
unset size
set terminal wxt