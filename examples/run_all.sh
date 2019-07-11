#!/bin/sh

top=`pwd`

for dir in `ls -F | grep /`
do
    echo ""
    echo "#######################################################################"
    echo "########  RUN EXAMPLE in $dir"
    echo "#######################################################################"
    echo ""

    cd $dir
    pwd

    if [ -f run.sh ]; then
        sh run.sh
    fi

    if [ -f check_results.sh ]; then
        sh check_results.sh
    fi

    cd $top
done
