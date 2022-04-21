
set style line 1 lt 1 lw 3 pt 2 ps 1.5
set style line 2 lt 2 lw 3 pt 2 ps 1.5
set style line 3 lt 3 lw 3 pt 4 ps 1.5
set style line 4 lt 4 lw 3 pt 4 ps 1.5


set key bottom
set zeroaxis lw 2

set xlabel "w_n"
set xrange[0:30]
set yrange[*:*]

set terminal postscript eps color enhanced "Times-Roman" 24


### Re
set ylabel "Re Sigma(iw_n)"
set output "sigma_re.eps"

plot "check/sigma.dat"u 1:2 title "ud 0" w lp ls 1\
, "" u 1:12 title "ud 1" w lp ls 2\
, "" u 1:22 title "ud 2" w lp ls 3\
, "" u 1:32 title "ud 3" w lp ls 4\


### Im
set ylabel "Im Sigma(iw_n)"
set output "sigma_im.eps"

plot "check/sigma.dat"u 1:3 title "ud 0" w lp ls 1\
, "" u 1:13 title "ud 1" w lp ls 2\
, "" u 1:23 title "ud 2" w lp ls 3\
, "" u 1:33 title "ud 3" w lp ls 4\


set output
