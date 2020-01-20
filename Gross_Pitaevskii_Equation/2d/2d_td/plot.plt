set terminal gif animate delay 10 optimize size 1600,1600
set output "data.gif"
set pm3d
# SETTINGS
dh          = 0.06
dt          = 0.01
N           = 100-1   # STEP COUNT
M           = N + 1
iter        = 5000 # ITERATION COUNT OF TIME
iter_output = 50    # SKIP COUNT IN THE ITERATION OF TIME
skip_output = 25    # SKIP COUNT OF SHOWING SPEED AND ESTIMATED TIME REMAINING
# OTHER VARIABLES USED BY SCRIPT
xmax        = real(N/2 + 0.5)*dh     # BOUNDARY OF X
data_num    = int(real(iter)/iter_output)
time_new    = 0.0
time_old    = 0.0
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
isdebug     = 1
stats "data.txt" u 3 nooutput
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
            splot "data_initial.txt" using 1:2:4 title "" with pm3d

            set zrange  [0:1]
            set cbrange [0:0.2]
            #set nosurface
            #set contour
            set title sprintf("Time development of Non-Linear Schroedinger Equation\n( T = %.3f )", dt*i*iter_output)
            splot "data.txt" using 1:2:3 every :::M*i::M*(i+1)-1 title "" with pm3d
            #unset contour
            #set surface

            set zrange [0:1]
            set cbrange [0:0.5]
            set title "Probability Profile of Initial Wave Function"
            splot "data_initial.txt" using 1:2:3 title "" with pm3d

            set zrange  [-50:50]
            set cbrange [-50:50]
            set title sprintf("External Potential")
            splot "data_potential.txt" using 1:2:4 every :::M*i::M*(i+1)-1 title "" with pm3d
        unset multiplot
    } else {
        if (STATS_max > 1) {
            set zrange  [0:STATS_max]
        } else {
            set zrange  [0:1]
        }
        set cbrange [0:STATS_max]
        #set nosurface
        #set contour
        set title sprintf("Time development of Non-Linear Schroedinger Equation\n( T = %.3f )", dt*i*iter_output)
        splot "data.txt" using 1:2:4 every :::M*i::M*(i+1)-1 title "" with pm3d
        #unset contour
        #set surface
    }
}


unset output
set terminal wxt