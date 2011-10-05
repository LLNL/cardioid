#!/bin/bash

# $Id$

searchPath_ddcMD="$DDCMD_SRC ../../ddcMD/src $HOME/src/ddcMD/src $HOME/ddcMD/src"

files_ddcMD="$@"

for dir in $searchPath_ddcMD
do
  if [ -r $dir/ddcMD.c ]
  then
    ddcPath=$dir
    break
  fi
done

if [ -z $ddcPath ]
then
    echo "ERROR: can't find ddcMD sources."
    echo "  search path:  $searchPath_ddcMD"
    exit 1
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
