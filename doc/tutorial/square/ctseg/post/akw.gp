set cbrange [0:0.8]
set size 0.95, 1.0
set xtics (\
  "G"  0.0, \
  "X"  0.5000000000000003, \
  "M"  1.0000000000000007, \
  "G"  1.7071067811865512, \
  "G"  1.8778174593052066, \
  "X"  2.377817459305207, \
  "M"  2.8778174593052075, \
  "G"  3.584924240491758 \
  )
set pm3d map
#set pm3d interpolate 5, 5
unset key
set ylabel "Energy"
set cblabel "A(k,w)"
splot "akw.dat" u 1:2:(abs($3))
pause -1

set term push
set term png
set out 'akw.png'
replot
set out
set term pop
