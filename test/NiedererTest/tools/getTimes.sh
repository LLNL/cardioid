#!/bin/bash

for balancer in koradi grid
do
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
    
    case $dx in
    0.1)
        xMin = -15
        yMin = -35
        zMin = -100
        xMax = 14
        yMax = 34
        zMax = 99
        ;;
    0.2)
        xMin = -8
        yMin = -18
        zMin = -50
        xMax = 7
        yMax = 17
        zMax = 49
        ;;
    0.5)
        xMin = -3
        yMin = -7
        zMin = -20
        xMax = 2
        yMax = 6
        zMax = 19
        ;;
    *)
        echo undefined dx
        exit 1
    esac

    P1 = 


    echo $balancer $dx $dt `grep "$target" $lastSnap/activationTime* | awk '{print $5}'`


    popd > /dev/null
done
done
done
