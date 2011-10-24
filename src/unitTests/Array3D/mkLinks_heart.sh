#!/bin/bash

# $Id$

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
