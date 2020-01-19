set terminal gif animate delay 10 optimize size 800,800
set size ratio 1
set output "rotation.gif"
unset pm3d
# SETTINGS
dh          = 0.06
dt          = 0.0001
N           = 100-1   # STEP COUNT
M           = N + 1
iter        = 5000 # ITERATION COUNT OF TIME
iter_output = 50    # SKIP COUNT IN THE ITERATION OF TIME
skip_output = 25    # SKIP COUNT OF SHOWING SPEED AND ETA
# OTHER VARIABLES USED BY SCRIPT
xmax        = real(N/2 + 0.5)*dh     # BOUNDARY OF X
data_num    = int(real(iter)/iter_output)
time_new    = 0.0
time_old    = 0.0
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set zrange [-1:1]
set xlabel "X (Unit:Xs)"
set ylabel "Y (Unit:Xs)"
set zlabel "rot(j)"
do for [i=0: data_num-1] {
    if (i%skip_output == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(skip_output) / interval
        time_old = time_new
        print sprintf("%5d / %5d    SPD : %5.2f lines/s   Estimated Time Remaining : %5.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }
    set title sprintf("Rotation of Probability Current\n( T = %.3f )", dt*i*iter_output)
    splot "data_rotation.txt" every :::2*M*i::2*M*(i+1)-1 title "" with vectors
}
unset output
set terminal wxt