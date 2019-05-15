set term svg enhanced size 1000, 750 fsize 32
set output "sigma.svg"

set xlabel "Energy"
beta = 50.0

set xr [0:2]
set yr [0:3]

set xlabel "(w_n)^{0.5}"
set ylabel "- Im Sigma(i w_n)"

plot \
"check/sigma.dat" u (($1)**0.5):($1 > 0 ? -$3 : 1/0) t "dxy" w p, \
"check/sigma.dat" u (($1)**0.5):($1 > 0 ? -$11 : 1/0) t "dyz" w p, \
"check/sigma.dat" u (($1)**0.5):($1 > 0 ? -$19 : 1/0) t "dzx" w p, \
"sigma-PRL101-166405.txt" u (((2*$0+1)*pi/beta)**0.5):($1) t "PRL 101, 166405 (2008)" lc 'r' lw 3 w l
 
#pause -1
