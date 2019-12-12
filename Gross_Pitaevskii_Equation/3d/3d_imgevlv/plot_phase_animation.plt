# GENERATES ANIMATION TO DISPLAY CROSS SECTION IMAGE AT EACH Z VALUE

set terminal gif animate delay 10 optimize size 800,800
set output "Phase_anim.gif"
# SETTINGS
xmax        = 5      # BOUNDARY OF X
N           = 30-1    # STEP COUNT
M           = N + 1
iter        = 2*(N+1) # ITERATION COUNT OF TIME
iter_output = 50      # SKIP COUNT IN THE ITERATION OF TIME
skip_output = 25      # SKIP COUNT OF SHOWING SPEED AND ETA
# OTHER VARIABLES USED BY SCRIPT
dh          = xmax / real(N/2 + 0.5)
dt          = 0.01*dh*dh
data_num    = iter
time_new    = 0.0
time_old    = 0.0
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set zrange [-xmax:xmax]
set xlabel "X (unit:Xs)"
set ylabel "Y (unit:Xs)"
set zlabel "Z (unit:Xs)"
set zlabel rotate by 90
set cbrange [-3.14:3.14]
do for [i=0: 2*N] {
    if (i%skip_output == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(skip_output) / interval
        time_old = time_new
        print sprintf("%5d / %5d    SPD : %5.2f lines/s   ETA : %5.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }
    set title "Result Wave Function"
    if (i > N) {
        j = 2*N-i
        splot "phase.txt" using 1:2:3:5 every :::M*j::(M*(j+1)-1) title "" with pm3d
    } else {
        splot "phase.txt" using 1:2:3:5 every :::M*i::(M*(i+1)-1) title "" with pm3d
    }
}
unset output
set terminal wxt