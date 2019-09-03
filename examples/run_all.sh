#!/bin/sh

topdir=`pwd`

for dir in `ls -F | grep /`
do
    cd $dir

    if [ -f run.sh ]; then
        echo "running example in $dir"
        sh run.sh 1>run.out 2>run.err
        echo " status: $?"
    fi

#    if [ -f check_results.sh ]; then
#        sh check_results.sh 1>check.out 2>check.err
#        echo " check: $?"
#    fi

    cd $topdir
done
