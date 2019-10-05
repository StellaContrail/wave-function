set terminal gif animate delay 1 optimize size 640,480
set title "Time dependance of Gaussian wave packet under Harmonic Oscillator Potential"
set output "data.gif"
set xrange [-15:15]
set yrange [-1:1]
set grid
do for [i=0: 630-1: 15] {
    plot "data.txt" every :::i::i using 1:2 with lines title "REAL","data.txt" every :::i::i using 1:3 with lines title "IMAG", "data.txt" every :::i::i using 1:4 with lines title "ABSL", 0.5*x*x*0.5 title "V(x)"
    #plot "data.txt" every :::i::i using 1:2 with lines title "REAL"
    #plot "data.txt" every :::i::i using 1:4 with lines title "PROB"
}
unset output
set terminal wxt