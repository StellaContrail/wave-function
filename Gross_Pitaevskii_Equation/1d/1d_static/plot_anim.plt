set terminal gif animate delay 10 optimize size 640,480
set output "data.gif"
set xrange [-10:10]
set grid
time_old = 0
data_num = 89
do for [i=0: data_num-1] {
    if (i%10 == 0) {
        time_new = time(0.0)
        interval = time_new - time_old
        speed    = real(10) / interval
        time_old = time_new
        print sprintf("%5d / %5d    SPD : %5.2f lines/s   ETA : %5.2f sec", i, data_num, speed, (data_num-i+1)/speed)
    }
    set title sprintf("i = %d", i)
    plot "data_td_potential.txt" every :::i::i with lines title "V(x)"
}
unset output
set terminal wxt