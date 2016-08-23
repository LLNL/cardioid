#!/bin/csh
if ($#argv != 1) then
  echo "Usage: $0 src_dirname"
  echo "concat files from one directory to one file with sorted data"
  exit 0
endif

set dirname = $1
set tmp1 = 'tmp1'
set tmp2 = 'tmp2'
set header = 'header'

set names = (anatomy activationTime coarsened_Ca)

foreach name ($names)
  echo $name ...
  set name0 = $name\#000000
  echo $name0 ...
  if ( -e $dirname/$name0 ) then
    set newdirname = {$dirname}_single
    echo 'New directory:' $newdirname
    if( ! -e $newdirname ) then
      mkdir $newdirname
    endif
    set newfile = $newdirname/$name\#000000
    rm -f $newfile
    echo 'remove' $newfile
    
    foreach file ( `ls -1 $dirname/$name*` )
      echo 'file: ' $file ...
      set nlines = `wc -l $file | awk '{print $1}'`
      set nheader = `grep FILEHEADER $file |wc -l`
      if ($nheader>0) then
        set n = `grep '=' $file | wc -l`
        set n = `echo $n+5 |bc`
        set nlines = `echo $nlines-$n |bc`
        tail -n $nlines $file >! $tmp2
        head -n $n $file >! $header
      else
        cp $file $tmp2
      endif
      if ( -e $newfile ) then
        cat $newfile $tmp2 >! $tmp1
      else
        cp $tmp2 $tmp1
      endif
      mv -f $tmp1 $newfile
    end
    echo 'sort...'
    sort $newfile >! $tmp1
    echo 'New file: ' $newfile
    cat $header $tmp1 >! $newfile
    rm -f $header
    rm -f $tmp1
    rm -f $tmp2
  endif
end
