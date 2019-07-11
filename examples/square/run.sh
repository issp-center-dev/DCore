#!/bin/sh

# ---------------------------------
# check path and environment variables

. ../tools/tools.sh

# check command
check_command dcore_pre
check_command dcore
check_command dcore_check
check_command dcore_post

${NUM_PROC:=1}  # set 1 if not defined

# ---------------------------------
# create and move into a directory

mkdir -p results
cd results

# ---------------------------------
# DCore

ini=../dmft_square.ini

echo "running dcore_pre..."
dcore_pre $ini

echo "running dcore..."
dcore --np $NUM_PROC $ini

echo "running dcore_check..."
dcore_check $ini

echo "running dcore_post..."
dcore_post --np $NUM_PROC $ini
