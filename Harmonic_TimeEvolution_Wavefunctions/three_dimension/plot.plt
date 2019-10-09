set terminal gif animate delay 1 optimize size 640,640
set pm3d map
set cbrange [0:0.25]
set size ratio 1 1
set output "data.gif"
set ylabel "y"
set xlabel "x"
do for [i=0: 23] {
    a = i * 40
    b = a + 39
    t = 0.02*i
    set title sprintf("t=%3.2f", t)
    splot "data.txt" every :::a::b using 1:2:5
}
unset output