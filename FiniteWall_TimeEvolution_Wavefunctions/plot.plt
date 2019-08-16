set terminal gif animate delay 5 optimize size 640,480
set output "data.gif"
set xrange [-40:40]
set yrange [-1:1]
set yzeroaxis
do for [i=0: 1050-1: 10] {
    #plot "data.txt" every :::i::i using 1:2 with lines title "REAL","data.txt" every :::i::i using 1:3 with lines title "IMAG", "data.txt" every :::i::i using 1:4 with lines title "PROB"
    plot "data.txt" every :::i::i using 1:4 with lines title "PROB"
}
unset output
set terminal wxt