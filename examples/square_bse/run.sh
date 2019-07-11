#!/bin/sh

# ---------------------------------
# check path and environment variables

which dcore
if [ $? = 1 ]; then
    echo "dcore not found" >&2
    exit 1
fi

which BSE
if [ $? = 1 ]; then
    echo "BSE not found" >&2
    exit 1
fi

which pomerol2dcore
if [ $? = 1 ]; then
    echo "pomerol2dcore not found" >&2
    exit 1
fi

if [ -z $BSE_DIR ]; then
    echo "Environment variable BSE_DIR not defined" >&2
    exit 1
fi

if [ -z $NUM_PROC ]; then
    echo "warning: NUM_PROC is not defined. Default value 1 will be used." >&2
fi

# ---------------------------------
# create and move into a directory

mkdir -p results
cp *in results
cd results

# ---------------------------------
# set variables

export PYTHONPATH=$BSE_DIR/python:$PYTHONPATH
echo "PYTHONPATH = $PYTHONPATH"

dir_script=`dirname $0`

ini=dmft_square.in
seedname=square

# ---------------------------------
# DCore

echo ""
echo "##################"
echo " DCore"
echo "##################"
echo ""

echo "running dcore_pre..."
dcore_pre $ini

echo "running dcore..."
dcore --np ${NUM_PROC:=1} $ini

echo "running dcore_check..."
dcore_check $ini

echo "running dcore_bse..."
dcore_bse --np ${NUM_PROC:=1} $ini

# ---------------------------------
# BSE

echo ""
echo "##################"
echo " BSE"
echo "##################"
echo ""

# print latest commit of git repository
$BSE_DIR/bin/misc/print_latest_commit.sh

# Generate q_path.dat
python $BSE_DIR/python/bse_tools/gen_qpath.py ${seedname}.h5 qpath.in

# Plot input to BSE
$BSE_DIR/python/plot/plot_bse_input.py

# BSE
python $BSE_DIR/python/bse_tools/bse_tool.py -s BSE -i dmft_bse.h5 -q q_path.dat
python $BSE_DIR/python/bse_tools/plot_chiq_path.py q_path.dat chiq_eigen.dat
python $BSE_DIR/python/bse_tools/plot_chiq_path.py q_path.dat chi0q_eigen.dat --mode='chi0'

# RPA
python $BSE_DIR/python/bse_tools/bse_tool.py -s BSE -i dmft_bse.h5 -q q_path.dat -t 'rpa'
python $BSE_DIR/python/bse_tools/plot_chiq_path.py q_path.dat chiq_rpa_eigen.dat --mode='rpa'

# RRPA
python $BSE_DIR/python/bse_tools/bse_tool.py -s BSE -i dmft_bse.h5 -q q_path.dat -t 'rrpa'
python $BSE_DIR/python/bse_tools/plot_chiq_path.py q_path.dat chiq_rrpa_eigen.dat --mode='rrpa'
