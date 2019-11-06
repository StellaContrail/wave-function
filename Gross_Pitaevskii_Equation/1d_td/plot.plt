set terminal gif animate delay 10 optimize size 640,480
set output "data.gif"
set xrange [-10:10]
set yrange [-1:1]
set grid
xmax        = 10
N           = 128
dh          = 2.0*xmax/N
dt          = 0.4*dh*dh
iter        = 50000
iter_output = 50
data_num    = int(iter/iter_output)
time_new    = 0.0
time_old    = 0.0
do for [i=0: data_num-1] {
    if (i%50 == 0) {
        interval = time_new - time_old
        speed    = interval == 0.0 ? 0.0 : 50.0 / interval
        time_old = time_new
        time_new = time(0.0)
        print sprintf("%d / %d    SPD : %.2f lines/s   ETA : %.2f sec", i, data_num, speed, speed==0.0?0.0:(data_num-i+1)/speed)
    }
    set title sprintf("Time development of Non-Linear Schroedinger Equation (t=%.2f s)", dt*i)
    #plot "data.txt" every :::i::i using 1:2 with lines title "REAL","data.txt" every :::i::i using 1:3 with lines title "IMAG", "data.txt" every :::i::i using 1:4 with lines title "ABSL", 0.5*x*x title "V(x)"
    #plot "data.txt" every :::i::i using 1:2 with lines title "REAL"
    plot "data.txt" every :::i::i using 1:4 with lines title "PROB", 0.5*x*x title "V(x)"
}
unset output
set terminal wxt