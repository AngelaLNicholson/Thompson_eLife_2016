#!/bin/bash

#Author: BZ

#usage: bash multiRunpipeline.sh prefsFile1.txt prefsFile2.txt ... prefsFileN.txt
#purpose: takes any # of pipeline2.0 prefs files, and runs the pipeline on each of them

i=0
num_simultaneous=8 #number to run at the same time. Don't overdo it.
for var in "$@"
do
    echo "$var"
    let 'i += 1'
    mod=$(expr $i % $num_simultaneous)
    bash runPipeline.sh $var &
    if [ $mod = 0 ]; then
        wait
    fi
done
wait