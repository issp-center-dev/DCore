set xl 'Energy'
set yl 'DOS'

set yr [0:]
plot 'dos.dat' w l

set term push
set term png
set out 'dos.png'
rep
set out
set term pop
