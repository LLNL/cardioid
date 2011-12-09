#!/bin/bash

for balancer in koradi grid
do
echo $balancer
for dx in 0.1 0.2 0.5
do
for dt in 0.005 0.01 0.05
do

    dirname=run/${balancer}_dt${dt}_dx${dx}
    pushd $dirname > /dev/null
    lastSnap=`ls -d snapshot.???????? |tail -n 1`
    if [ ! -r $lastSnap/activationTime#000000 ]
    then
        popd > /dev/null
        continue
    fi
    #echo $dirname/$lastSnap
    
    echo -n "$dx $dt & "
    cat $lastSnap/activationTime#* | ../../tools/getTimes.py
    echo \\hline

    popd > /dev/null
done
done
done
