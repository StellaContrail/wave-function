set title "Tunneling Effect"
set ylabel "Probability"
set xlabel "X / x0"
set yrange [0:1]
set xrange [-5:5]
plot for [i=1:8] 'data.txt' every :::4*i::4*i using 1:4 w l title sprintf("WF %d", i)
replot "data_potential.txt" every :::0::0 w l title "POTENTIAL"
