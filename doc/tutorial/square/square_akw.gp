#set term png
#set output "akw.png"
set size 0.95, 1.0
set xtics (\
  "G"  0.0, \
  "X"  0.5, \
  "M"  1.0, \
  "G"  1.70710678119, \
  "G"  1.87781745931, \
  "X"  2.37781745931, \
  "M"  2.87781745931, \
  "G"  3.58492424049 \
  )
set pm3d map
#set pm3d interpolate 5, 5
unset key
set ylabel "Energy"
set cblabel "A(k,w)"
splot "square_akw.dat", \
"square_akw0.dat" u 1:($2-0.0):(0) every 5 w p lc 5
#pause -1
