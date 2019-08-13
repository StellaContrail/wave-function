set terminal gif animate delay 5 optimize size 640,480
set output "data.gif"
set xrange [-20:20]
set yrange [-1:1]
set grid
do for [i=0: 10000-1: 15] {
    plot "data.txt" every :::i::i using 1:2 with lines title "REAL","data.txt" every :::i::i using 1:3 with lines title "IMAG", "data.txt" every :::i::i using 1:4 with lines title "PROB"
    #plot "data.txt" every :::i::i using 1:2 with lines title "REAL"
    #plot "data.txt" every :::i::i using 1:4 with lines title "PROB"
}
unset output
set terminal wxt