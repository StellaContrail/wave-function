set terminal gif animate delay 10 optimize size 1600,800
set output "data_flux.gif"
set yrange [-10:10]
set xrange [-10:10]
# SETTINGS
xmax        = 10     # BOUNDARY OF X
N           = 50-1   # STEP COUNT
iter        = 10000 # ITERATION COUNT OF TIME
iter_output = 50    # SKIP COUNT IN THE ITERATION OF TIME
skip_output = 25    # SKIP COUNT OF SHOWING SPEED AND ETA
# OTHER VARIABLES USED BY SCRIPT
dh          = xmax / real(N/2 + 0.5)
dt          = 0.01*dh*dh
data_num    = int(real(iter)/iter_output)
time_new    = 0.0
time_old    = 0.0
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set xlabel "X (Unit:Xs)"
set ylabel "Y (Unit:Xs)"
set zlabel rotate by 90
do for [i=0: data_num-1] {
    if (i%skip_output == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(skip_output) / interval
        time_old = time_new
        print sprintf("%5d / %5d    SPD : %5.2f lines/s   ETA : %5.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }
    set multiplot layout 1,2
        set zlabel "Probability"
        set pm3d
        set zrange [0:1]
        set cbrange [0:0.2]
        set title sprintf("Time development of Wave Function\n( T = %.3f )", dt*i)
        splot "data_projection.txt" using 1:2:3 every :::50*i::50*(i+1)-1 title "" with pm3d

        set zlabel "Z (Unit:Xs)"
        set zrange [-15:15]
        unset pm3d
        set title "Probability Current"
        splot "data_flux.txt" every :::2500*i::(2500*(i+1)-1) title "" with vectors
    unset multiplot
}
unset output
set terminal wxt