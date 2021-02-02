set size 0.95, 1.0
set xtics (\
  "G"  0.0, \
  "X"  0.5000000000000003 \
  )
set pm3d map
#set pm3d interpolate 5, 5
unset key
set ylabel "Energy"
set cblabel "A(k,w)"
splot "test_akw.dat"
pause -1
