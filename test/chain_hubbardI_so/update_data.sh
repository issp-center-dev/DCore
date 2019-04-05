rm -rf ref
mkdir ref

dcore_pre dmft.ini
cp test.h5 ref/

dcore --np 1 dmft.ini
cp test.out.h5 ref/

dcore_post --np 1 dmft.ini
cp test*.dat ref/
