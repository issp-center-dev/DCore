set term png
set output "results/test_akw.png"
set size 0.95, 1.0
set xtics (\
  "G"  0.0, \
  "X"  0.5000000000000003, \
  "G"  0.5500000000000004, \
  "X"  1.0500000000000007 \
  )
set pm3d map
#set pm3d interpolate 5, 5
unset key
set ylabel "Energy"
set cblabel "A(k,w)"
splot "results/post/test_akw.dat" u 1:2:(abs($3))
