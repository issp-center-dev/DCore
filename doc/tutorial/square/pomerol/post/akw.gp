set terminal unknown
load "square_akw.gp"


#set terminal postscript eps color enhanced "Times-Roman" 24
#set output "akw.eps"
#set outpu "| epstopdf -f -o=akw.pdf"

set terminal png
set output "akw.png"

#set palette rgbformulae 21, 22, 23  # hot
#set palette defined (0 "white", 1 "blue", 2 "orange")  # white background

set cbrange[0:0.8]
set zrange[*:*]

replot
set output
