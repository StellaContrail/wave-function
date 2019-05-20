set terminal gif animate delay 5 optimize size 640,480
set output "data.gif"
do for [i=0: 500] {
    plot [-20:20][-1:1] "data.txt" every :::i::i using 2:3 with lines
}
unset output