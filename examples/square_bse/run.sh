#!/bin/sh

# ---------------------------------
# check path and environment variables

. ../tools/tools.sh

# check command
check_command dcore
check_command BSE
check_command pomerol2dcore

# check environment variable
check_var BSE_DIR

${NUM_PROC:=1}  # set 1 if not defined

# ---------------------------------
# create and move into a directory

mkdir -p results
#cp *in results
cd results

# ---------------------------------
# set variables

export PYTHONPATH=$BSE_DIR/python:$PYTHONPATH
echo "PYTHONPATH = $PYTHONPATH"

dir_script=`dirname $0`

ini=../dmft_square.in
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
dcore --np $NUM_PROC $ini

echo "running dcore_check..."
dcore_check $ini

echo "running dcore_bse..."
dcore_bse --np $NUM_PROC $ini

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
python $BSE_DIR/python/bse_tools/gen_qpath.py ${seedname}.h5 ../qpath.in

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
