set size 0.95, 1.0
set xtics (\
  "{/Symbol G}"  0.0, \
  "X"  0.5, \
  "M"  1.0, \
  "{/Symbol G}"  1.70710678119, \
  "{/Symbol G}"  1.87781745931, \
  "X"  2.37781745931, \
  "M"  2.87781745931, \
  "{/Symbol G}"  3.58492424049 \
  )
set pm3d map
#set pm3d interpolate 5, 5
unset key
set ylabel "Energy (eV)"
#set cblabel "A(k,w)"


set terminal postscript eps color enhanced "Times-Roman" 24
#set output "square_akw_paper.eps"
set outpu "| epstopdf -f -o=square_akw_paper.pdf"

#set terminal png size  1600, 1200
#set output "akw.png"

#set palette rgbformulae 21, 22, 23  # hot
#set palette defined (0 "white", 1 "blue", 2 "orange")  # white background

set cbrange[0:0.8]
set zrange[*:*]



splot "square_akw.dat"
set output
