#!/bin/bash

USER_ID=${LOCAL_UID:-9001}
GROUP_ID=${LOCAL_GID:-9001}

echo "Starting with UID : $USER_ID, GID: $GROUP_ID"
ls -al /home/
useradd -u $USER_ID -o -m user -s /bin/bash
groupmod -g $GROUP_ID user
export HOME=/home/user

cd /var/dcoretest
#git clean -xdf

pip3 install pytest matplotlib
pip3 install .[dev]

chown -R user /var/dcoretest
chgrp -R user /var/dcoretest
#gosu user bash -c 'source /opt/triqs/share/triqsvars.sh;MPIRUN=mpirun NUM_PROC=2 ls examples'
#gosu user bash -c 'source /opt/triqs/share/triqsvars.sh;cd examples;MPIRUN=mpirun NUM_PROC=2 DCORE_CHECK_DEFAULT_EXT=png sh run_ci.sh'
#gosu user bash -c '/bin/bash'