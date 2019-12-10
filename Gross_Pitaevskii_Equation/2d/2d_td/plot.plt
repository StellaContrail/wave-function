set terminal gif animate delay 10 optimize size 1600,1600
set output "data.gif"
set pm3d
# SETTINGS
xmax        = 5     # BOUNDARY OF X
N           = 30-1   # STEP COUNT
M           = N + 1
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
do for [i=0: data_num-1] {
    if (i%skip_output == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(skip_output) / interval
        time_old = time_new
        print sprintf("%5d / %5d    SPD : %5.2f lines/s   ETA : %5.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }
    set multiplot layout 2,2 scale 1,1
        set zrange [-1:1]
        set cbrange [-1:1]
        set title "Real Part Profile of Initial Wave Function"
        splot "data_initial.txt" using 1:2:4 title "" with pm3d

        set zrange  [0:1]
        set cbrange [0:1]
        set title sprintf("Time development of Non-Linear Schroedinger Equation\n( T = %.3f )", dt*i*iter_output)
        splot "data.txt" using 1:2:3 every :::M*i::M*(i+1)-1 title "" with pm3d


        set zrange [0:1]
        set cbrange [0:1]
        set title "Probability Profile of Initial Wave Function"
        splot "data_initial.txt" using 1:2:3 title "" with pm3d

        set zrange  [-5:5]
        set cbrange [-5:5]
        set title sprintf("External Potential")
        splot "data_potential.txt" using 1:2:3 every :::M*i::M*(i+1)-1 title "" with pm3d
    unset multiplot
}
unset output
set terminal wxt