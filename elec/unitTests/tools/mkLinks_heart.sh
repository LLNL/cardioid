#!/bin/bash

# $Id: mkLinks_heart.sh 32 2011-10-24 22:47:06Z richards12 $

searchPath="../../"

files="$@"

for file in $files
do
    echo $file
    if [ ! -r $file ]
    then
        ln -s $searchPath/$file
    fi
done
