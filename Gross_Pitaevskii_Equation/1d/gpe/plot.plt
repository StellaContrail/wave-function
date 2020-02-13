set terminal gif animate delay 10 optimize size 800,600
set output "data.gif"
load "gnuplot_vars.txt"
# SETTINGS
M           = N + 1
skip_output = 25    # SKIP COUNT OF SHOWING SPEED AND ESTIMATED TIME REMAINING
data_num    = int(real(iter)/iter_output)
time_new    = 0.0
time_old    = 0.0
stats fn_wavefunction_real_result u 3 nooutput
set xrange [-xmax:xmax]
set xrange [-10:10]
set yrange [0:0.2]
set xlabel "X/a0"
set ylabel "Density"
set grid
isdebug     = 1
set xtics 5
set ytics 0.05

set terminal png size 800, 600
set output "data_0.png"
set title sprintf("Density\n( T = %.3f )", dt*0*iter_output)
plot fn_wavefunction_real_result using 1:2 every :::0::0 title "" with lines

set output "data_1.png"
set title sprintf("Density\n( T = %.3f )", dt*((data_num-1)/2)*iter_output)
plot fn_wavefunction_real_result using 1:2 every :::(data_num-1)/2::(data_num-1)/2 title "" with lines

set output "data_2.png"
set title sprintf("Density\n( T = %.3f )", dt*(data_num-1)*iter_output)
plot fn_wavefunction_real_result using 1:2 every :::data_num-1::data_num-1 title "" with lines
pause -1 'test'

do for [i=0: data_num-1] {
    if (i%skip_output == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(skip_output) / interval
        time_old = time_new
        print sprintf("%5d / %5d    SPD : %5.2f lines/s   Estimated Time Remaining : %5.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }

    if (isdebug == 0) {
        set multiplot layout 2,2 scale 1,1
            set zrange [-1:1]
            set cbrange [-1:1]
            set title "Real Part Profile of Initial Wave Function"
            splot fn_wavefunction_imaginary_result using 1:2:4 title "" with pm3d

            set zrange  [0:1]
            set cbrange [0:1]
            #set nosurface
            #set contour
            set title sprintf("Time development of Non-Linear Schroedinger Equation\n( T = %.3f )", dt*i*iter_output)
            splot fn_wavefunction_real_result using 1:2:3 every :::M*i::M*(i+1)-1 title "" with pm3d
            #unset contour
            #set surface

            set zrange [0:1]
            set cbrange [0:1]
            set title "Probability Profile of Initial Wave Function"
            splot fn_wavefunction_imaginary_result using 1:2:3 title "" with pm3d

            unset zrange
            unset cbrange
            set title sprintf("External Potential")
            splot fn_potential_imaginary using 1:2:4 every :::M*i::M*(i+1)-1 title "" with pm3d
        unset multiplot
    } else {
        set title sprintf("Density\n( T = %.3f )", dt*i*iter_output)
        plot fn_wavefunction_real_result using 1:2 every :::i::i title "" with lines#, fn_potential_imaginary using 1:2 title "" with lines
    }
}
unset output
set terminal wxt
