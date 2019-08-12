set terminal gif animate delay 10 optimize size 640,480
set output "data.gif"
set xrange [-25:25]
set yrange [-1:1]
do for [i=0: 1000-1: 1] {
    plot "data.txt" every :::i::i using 1:2 with lines title "REAL","data.txt" every :::i::i using 1:3 with lines title "IMAG", "data.txt" every :::i::i using 1:4 with lines title "PROB"
}
unset output