set term png
set output "dos.png"

set xlabel "Energy"
set ylabel "DOS"
plot "square_dos.dat" w l
