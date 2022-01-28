#!/bin/sh

# ---------------------------------
# check path and environment variables

. ../tools/tools.sh

# check command
check_command dcore_pre
check_command dcore
check_command dcore_check
check_command dcore_post

# set NUM_PROC=1 if not defined
set_num_proc

# ---------------------------------
# for spin_orbit=False
# ---------------------------------
# create and move into a directory

rm -rf results
mkdir results
cp square_2orb_hr.dat results/
cd results || exit

# ---------------------------------
# DCore

ini=../square_2orb.ini

echo "running dcore_pre..."
dcore_pre $ini || exit

echo "running dcore..."
dcore --np $NUM_PROC $ini || exit

echo "running dcore_check..."
dcore_check $ini || exit

echo "running dcore_post..."
dcore_post --np $NUM_PROC $ini || exit

cd ..

# ---------------------------------
# for spin_orbit=True
# ---------------------------------
# create and move into a directory

rm -rf results_so
mkdir results_so
cp square_2orb_so_hr.dat results_so/
cd results_so || exit

# ---------------------------------
# DCore

ini=../square_2orb_so.ini

echo "running dcore_pre..."
dcore_pre $ini || exit

echo "running dcore..."
dcore --np $NUM_PROC $ini || exit

echo "running dcore_check..."
dcore_check $ini || exit

echo "running dcore_post..."
dcore_post --np $NUM_PROC $ini || exit
