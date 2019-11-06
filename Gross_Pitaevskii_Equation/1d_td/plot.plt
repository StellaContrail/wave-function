set terminal gif animate delay 1 optimize size 640,480
set output "data.gif"
set xrange [-10:10]
set yrange [0:1]
set grid
do for [i=0: 100-1] {
    plot "data.txt" every :::i::i using 1:2 with lines title "ABSL"
}
unset output
set terminal wxt