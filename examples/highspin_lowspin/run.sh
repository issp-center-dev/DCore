#!/bin/sh

# ---------------------------------
# check path and environment variables

. ../tools/tools.sh

# check command
check_command dcore
check_command dcore_bse
check_command BSE
check_command pomerol2dcore

# check environment variable
check_var BSE_DIR

# set NUM_PROC=1 if not defined
set_num_proc

# ---------------------------------
# create and move into a directory

rm -rf results
mkdir results
cd results
cp ../eigenvec.in .

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
check_status

echo "running dcore..."
dcore --np $NUM_PROC $ini
check_status

echo "running dcore_check..."
dcore_check $ini
check_status

echo "running dcore_bse..."
dcore_bse --np $NUM_PROC $ini
check_status

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
python3 $BSE_DIR/python/bse_tools/gen_qpath.py ${seedname}.h5 ../qpath.in
check_status

# BSE
mpirun -np $NUM_PROC python3 $BSE_DIR/python/bse_tools/bse_tool.py --toml ../bse.toml
check_status
mpirun -np $NUM_PROC python3 $BSE_DIR/python/bse_tools/bse_post.py --toml ../bse.toml
check_status

# Plot BSE results
python3 $BSE_DIR/python/bse_tools/plot_chiq_path.py q_path.dat chi_q_eigen.dat
python3 $BSE_DIR/python/bse_tools/plot_chiq_path.py q_path.dat chi0_q_eigen.dat --mode='chi0'
