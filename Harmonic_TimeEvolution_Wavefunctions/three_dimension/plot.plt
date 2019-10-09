set terminal gif animate delay 1 optimize size 640,640
set pm3d map
set cbrange [0:0.25]
set xrange [-15.2:16]
set yrange [-15.2:16]
set size ratio 1 1
set output "data.gif"
set ylabel "Y"
set xlabel "X"
set palette cubehelix start 0.5 cycles -1.5 saturation 1
m = 1200
do for [i=0: m/25-1] {
    a = i * 40
    b = a + 39
    t = 0.02*(i*25)
    set title sprintf("T = %3.2f", t)
    splot "data.txt" every :::a::b using 1:2:3
}
unset output