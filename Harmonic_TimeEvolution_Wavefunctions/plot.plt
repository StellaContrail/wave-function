set terminal gif animate delay 1 optimize size 640,480
set output "data.gif"
set xrange [-20:20]
set yrange [-1:1]
do for [i=0: 2500: 100] {
    plot "data.txt" every :::i::i using 2:3 with lines title "REAL","data.txt" every :::i::i using 2:4 with lines title "IMAG","data.txt" every :::i::i using 2:5 with lines title "PROB"
}
unset output