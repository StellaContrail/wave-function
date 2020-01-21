set terminal gif animate delay 10 optimize size 800,800
set size ratio 1
set output "current.gif"
# SETTINGS
dh          = 0.16
dt          = 0.002
N           = 100-1   # STEP COUNT
M           = N + 1
iter        = 1000 # ITERATION COUNT OF TIME
iter_output = 10    # SKIP COUNT IN THE ITERATION OF TIME
skip_output = 25    # SKIP COUNT OF SHOWING SPEED AND ETA
# OTHER VARIABLES USED BY SCRIPT
xmax        = real(N/2 + 0.5)*dh     # BOUNDARY OF X
data_num    = int(real(iter)/iter_output)
time_new    = 0.0
time_old    = 0.0
set xrange [-xmax:xmax]
set yrange [-xmax:xmax]
set xlabel "X (Unit:Xs)"
set ylabel "Y (Unit:Xs)"
set view 0, 0, 1, 1
set palette rgbformulae 22,13,10
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
    splot "data.txt" using 1:2:3 every :::M*i::M*(i+1)-1 title "", "data_flux.txt" using 1:2:(0):($3*15):($4*15):(0) every :::M*i::M*(i+1)-1 title "" with vectors linecolor rgb "#000000"
    unset pm3d
    
    #set title sprintf("Probability current\n( T = %.3f )", dt*i*iter_output)
}
unset output
set terminal wxt
 set view 60, 30, 1, 1