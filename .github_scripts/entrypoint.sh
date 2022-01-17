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
pip3 install git+https://github.com/shinaoka/dcorelib.git
pip3 install .

chown -R user /var/dcoretest
chgrp -R user /var/dcoretest
gosu user bash -c 'source /opt/triqs/share/triqsvars.sh;pytest tests/non-mpi/*/*.py'
gosu user bash -c 'source /opt/triqs/share/triqsvars.sh;pytest tests/non-mpi/test*.py'
gosu user bash -c 'source /opt/triqs/share/triqsvars.sh;mpirun -np 2 pytest tests/mpi/*/*.py'
echo "TEST DONE"
#echo "Use Ctrl+C to stop the container!"
