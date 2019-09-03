#!/bin/sh

# ---------------------------------
# check path and environment variables

. ../tools/tools.sh

# check command
check_command dcore_pre
check_command dcore
check_command dcore_check
check_command hybmat

# check environment variable
check_var MPIRUN

# set NUM_PROC=1 if not defined
set_num_proc

# ---------------------------------
# create and move into a directory

rm -rf results
mkdir results
cd results

# ---------------------------------
# DCore

ini=../dmft_bethe.ini

echo "running dcore_pre..."
dcore_pre $ini
check_status

echo "running dcore..."
dcore --np $NUM_PROC $ini
check_status

echo "running dcore_check..."
dcore_check $ini
check_status
