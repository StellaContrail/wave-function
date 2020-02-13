set terminal gif animate delay 10 optimize size 800,800
set size ratio 1
set output "current.gif"
load "gnuplot_vars.txt"
# SETTINGS
M           = N + 1
skip_output = 25    # SKIP COUNT OF SHOWING SPEED AND ETA
# OTHER VARIABLES USED BY SCRIPT
data_num    = int(real(iter)/iter_output)
time_new    = 0.0
time_old    = 0.0
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set xlabel "X (Unit:Xs)"
set ylabel "Y (Unit:Xs)"
set view 0, 0, 1, 1
set palette rgbformulae 22,13,10
stats fn_wavefunction_real_result using 3 nooutput
set zrange [0:1]
set cbrange [0:STATS_max]
stats fn_current_real_result using 3:4 nooutput
SCALE = 1.0 / sqrt(STATS_max_x**2 + STATS_max_y**2)
do for [i=0: data_num-1] {
    if (i%skip_output == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(skip_output) / interval
        time_old = time_new
        print sprintf("%5d / %5d    SPD : %5.2f lines/s   Estimated Time Remaining : %5.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }
    set pm3d map
    set title sprintf("Probability current\n( T = %.3f )", dt*i*iter_output)
    splot fn_wavefunction_real_result using 1:2:3 every :::M*i::M*(i+1)-1 title "", fn_current_real_result using 1:2:(0):($3*SCALE):($4*SCALE):(0) every :::M*i::M*(i+1)-1 title "" with vectors linecolor rgb "#000000"
    unset pm3d
    
    #set title sprintf("Probability current\n( T = %.3f )", dt*i*iter_output)
}
unset output
set terminal wxt
 set view 60, 30, 1, 1