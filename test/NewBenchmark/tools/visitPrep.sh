#!/bin/bash

if [ $# -ne 1 ]
then
    echo "Must supply filename"
    exit -1
fi


for dir in snapshot.????????
do
    if [ -e $dir.ddcMD ]
    then
        continue
    fi

    echo "ddcMD FILEHEADER {files=$1;}" > $dir.ddcMD
done
