#!/bin/sh

# ---------------------------------
# check path and environment variables

. ../tools/tools.sh

# check command
check_command dcore
check_command dcore_bse
check_command pomerol2dcore
check_command BSE
check_command bse_tool.py
check_command bse_post.py
check_command gen_qpath.py
check_command plot_chiq_path.py

# set NUM_PROC=1 if not defined
set_num_proc

# ---------------------------------
# create and move into a directory

rm -rf results
mkdir results
cd results || exit

# ---------------------------------
# DCore

ini=../dmft_square.in

echo ""
echo "##################"
echo " DCore"
echo "##################"
echo ""

echo "running dcore_pre..."
dcore_pre $ini || exit

echo "running dcore..."
dcore --np $NUM_PROC $ini || exit

echo "running dcore_check..."
dcore_check $ini || exit

echo "running dcore_bse..."
dcore_bse --np $NUM_PROC $ini || exit

# ---------------------------------
# BSE

echo ""
echo "##################"
echo " BSE"
echo "##################"
echo ""

# Print version
bse_tool.py --version || exit

# Generate q_path.dat
gen_qpath.py $ini ../qpath.in || exit

# BSE
mpirun -np $NUM_PROC bse_tool.py ../bse.in || exit
bse_post.py ../bse.in || exit

# Plot BSE results
plot_chiq_path.py q_path.dat chi_q_eigen.dat
plot_chiq_path.py q_path.dat chi0_q_eigen.dat --mode='chi0'
plot_chiq_path.py q_path.dat chi_q_rpa_eigen.dat --mode='rpa'
plot_chiq_path.py q_path.dat chi_q_rrpa_eigen.dat --mode='rrpa'
