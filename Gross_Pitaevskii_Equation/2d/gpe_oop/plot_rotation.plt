set terminal gif animate delay 10 optimize size 800,800
set size ratio 1
set output "rotation.gif"
set pm3d map
set view 0, 0, 1, 1
load "gnuplot_vars.txt"
stats fn_rotation_real_result u 3 nooutput
# SETTINGS
M           = N + 1
skip_output = 25    # SKIP COUNT OF SHOWING SPEED AND ETA
# OTHER VARIABLES USED BY SCRIPT
data_num    = int(real(iter)/iter_output)
time_new    = 0.0
time_old    = 0.0
set grid
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set zrange [-10:10]
set cbrange [-10:10]
set xlabel "X (Unit:Xs)"
set ylabel "Y (Unit:Xs)"
set zlabel "rot(j)"
set zlabel rotate by 90
do for [i=0: data_num-1] {
    if (i%skip_output == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(skip_output) / interval
        time_old = time_new
        print sprintf("%5d / %5d    SPD : %5.2f lines/s   Estimated Time Remaining : %5.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }
    set title sprintf("Rotation of Probability Current\n( T = %.3f )", dt*i*iter_output)
    splot fn_rotation_real_result every :::2*M*i::2*M*(i+1)-1 u 1:2:4 title "" with pm3d
}
unset output
set terminal wxt