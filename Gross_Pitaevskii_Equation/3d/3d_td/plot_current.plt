set terminal gif animate delay 10 optimize size 640,480
set output "data_current.gif"
set yrange [-1:1]
set grid
# SETTINGS
xmax        = 16     # BOUNDARY OF X
N           = 257   # STEP COUNT
iter        = 50000 # ITERATION COUNT OF TIME
iter_output = 50    # SKIP COUNT IN THE ITERATION OF TIME
skip_output = 25    # SKIP COUNT OF SHOWING SPEED AND ETA
# OTHER VARIABLES USED BY SCRIPT
dh          = real(2*xmax)/N
dt          = 0.4*dh*dh
data_num    = int(real(iter)/iter_output)
time_new    = 0.0
time_old    = 0.0
set xrange [-xmax:xmax]
do for [i=0: data_num-1] {
    if (i%skip_output == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(skip_output) / interval
        time_old = time_new
        print sprintf("%d / %d    SPD : %.2f lines/s   ETA : %.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }
    # [NOTE]
    # When "GD Warning: one parameter to a memory allocation multiplication is negative or zero, failing operation gracefully" appears,
    # it is a sign that a result of some multiplication is not shown properly (or change of degit isn't shown entirely)
    set title sprintf("Time development of Probability current of 1D TDGPE\n( T = %.3f s )", dt*i)
    #plot "data.txt" every :::i::i using 1:2 with lines title "REAL","data.txt" every :::i::i using 1:3 with lines title "IMAG", "data.txt" every :::i::i using 1:4 with lines title "ABSL", 0.5*x*x title "V(x)"
    #plot "data.txt" every :::i::i using 1:2 with lines title "REAL"
    plot "data_current.txt" every :::i::i with lines title "CURRENT"
}
unset output
set terminal wxt