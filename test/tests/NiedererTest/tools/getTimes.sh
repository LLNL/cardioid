#!/bin/bash

for balancer in koradi grid
do
echo $balancer
for dx in 0.10 0.20 0.50
do
for dt in 0.005 0.010 0.050
do

    dirname=run/${balancer}_dt${dt}_dx${dx}
    pushd $dirname > /dev/null
    lastSnap=`ls -d snapshot.???????????0 |tail -n 1`
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
