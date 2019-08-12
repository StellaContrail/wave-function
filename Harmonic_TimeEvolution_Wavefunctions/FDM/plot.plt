set terminal gif animate delay 1 optimize size 640,480
set output "data.gif"
set xrange [-25:25]
set yrange [-1:1]
do for [i=0: 2000-1: 100] {
    plot "data.txt" every :::i::i using 1:2 with lines title "REAL","data.txt" every :::i::i using 1:3 with lines title "IMAG"
}
unset output