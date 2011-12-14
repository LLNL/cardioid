#!/bin/bash

# $Id$

searchPath_ddcMD="$DDCMD_SRC $HOME/src/ddcMD/src $HOME/ddcMD/src ddcMD_files/src"

files_ddcMD="$@"

for dir in $searchPath_ddcMD
do
  if [ -r $dir/object.c ]
  then
    ddcPath=$dir
    break
  fi
done

if [ -z $ddcPath ]
then
    if [ -e ddcMD_files_r????.tgz ]
    then
        tar -xzf ddcMD_files_r????.tgz 
        ddcPath=ddcMD_files/src
    else
        echo "ERROR: can't find ddcMD sources."
        echo "  search path:  $searchPath_ddcMD"
        exit 1
    fi
else
    echo "Using ddcMD sources from $ddcPath"
fi

for file in $files_ddcMD
do
    echo $file
    if [ ! -r $file ]
    then
	ln -s $ddcPath/$file
    fi
done
