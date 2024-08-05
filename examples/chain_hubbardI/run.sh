#!/bin/sh

# ---------------------------------
# check path and environment variables

. ../tools/tools.sh

# check command
check_command dcore_pre
check_command dcore
check_command dcore_check
check_command dcore_anacont
check_command dcore_spectrum

# set NUM_PROC=1 if not defined
set_num_proc

# ---------------------------------
# create and move into a directory

rm -rf results
mkdir results
cd results || exit

# ---------------------------------
# DCore

ini=../dmft.ini

echo "running dcore_pre..."
dcore_pre $ini || exit

echo "running dcore..."
dcore --np $NUM_PROC $ini || exit

echo "running dcore_check..."
dcore_check $ini || exit

echo "running dcore_anacont..."
dcore_anacont $ini || exit

echo "running dcore_spectrum..."
dcore_spectrum --np $NUM_PROC $ini || exit

echo "plotting akw..."
cd ..
gnuplot test_akw.gp
