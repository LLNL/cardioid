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
        target="  14    34    99  "
        ;;
    0.2)
        target="  7    17    49 "
        ;;
    0.5)
        target="   2     6    19 "
        ;;
    *)
        echo undefined dx
        exit 1
    esac

    echo $balancer $dx $dt `grep "$target" $lastSnap/activationTime* | awk '{print $5}'`


    popd > /dev/null
done
done
done
