#set terminal postscript enhanced color 28
set terminal png 24

set output 'fit.png'
set log y
plot '<grep epoch work/sparse_fit-D50/output-wb0.txt' u 2:($6) t 'root sequared error' w lp, \
'<grep epoch work/sparse_fit-D50/output-wb0.txt' u 2:($8/$6) t 'relative diff of root sequared error' w lp
