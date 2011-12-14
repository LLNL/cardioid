#!/bin/bash

# $Id$

searchPath_ddcMD="$DDCMD_SRC $HOME/src/ddcMD/src $HOME/ddcMD/src"

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
    echo "ERROR: can't find ddcMD sources."
    echo "  search path:  $searchPath_ddcMD"
    exit 1
else
    echo "Using ddcMD sources from $ddcPath"
fi

rm -rf ddcMD_files
mkdir -p ddcMD_files/src

for file in $files_ddcMD
do
    cp $ddcPath/$file ddcMD_files/src
done

svnVersion=`(cd $ddcPath && svnversion)`

tar -czf ddcMD_files_r${svnVersion}.tgz ddcMD_files
