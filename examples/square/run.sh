#!/bin/sh

# ---------------------------------
# check path and environment variables

which dcore
if [ $? = 1 ]; then
    echo "dcore not found" >&2
    exit 1
fi

if [ -z $NUM_PROC ]; then
    echo "warning: NUM_PROC is not defined. Default value 1 will be used." >&2
fi

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
dcore --np ${NUM_PROC:=1} $ini

echo "running dcore_check..."
dcore_check $ini

echo "running dcore_post..."
dcore_post --np ${NUM_PROC:=1} $ini
