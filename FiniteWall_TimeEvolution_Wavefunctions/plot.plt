set terminal gif animate delay 5 optimize size 640,480
set output "data.gif"
set title "Tunneling Effect by traveling Gaussian Wave packet"
set xrange [-50:50]
set yrange [0:0.5]
set yzeroaxis
do for [i=0: 1000-1: 25] {
    #plot "data.txt" every :::i::i using 1:2 with lines title "REAL","data.txt" every :::i::i using 1:3 with lines title "IMAG", "data.txt" every :::i::i using 1:4 with lines title "PROB"
    plot "data.txt" every :::i::i using 1:4 with lines title "PROB"
}
unset output
set terminal wxt