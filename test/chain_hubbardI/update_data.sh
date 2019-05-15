rm -rf ref
mkdir ref

dcore_pre dmft.ini
dcore --np 1 dmft.ini
dcore_post --np 1 dmft.ini

cp test.h5 ref/
cp test.out.h5 ref/
cp test*.dat ref/
