#!/bin/bash

# To Do

# 1. Might want to add a loop over different numbers of tasks to check
# that the results don't change with the decomposition.  

exe=
pool=
maxTime=
bank=
nNodes=
nTasks=


for dt in 0.05 0.01 0.005
do
for dx in 0.5 0.2 0.1
do
for balancer in koradi grid
do

  dirname=run/${balancer}_dt${dt}_dx${dx}
  if [ -d $dirname ] 
  then
      continue
  fi

  echo making $dirname
  mkdir -p $dirname

  case $dx in
      0.5)
      stimBox=4
      ;;
      0.2)
      stimBox=8
      ;;
      0.1)
      stimBox=16
      ;;
      *)
      echo ERROR: undefined dx
      exit -1
  esac

  case $dt in
      0.05)
      maxLoop=20000
      ;;
      0.01)
      maxLoop=100000
      ;;
      0.005)
      maxLoop=200000
      ;;
      *)
      echo ERROR: undefined dt
      exit -1
  esac


  cat tools/object.data.proto \
      | sed s/XX_DT_XX/$dt/ \
      | sed s/XX_DX_XX/$dx/ \
      | sed s/XX_MAXLOOP_XX/$maxLoop/ \
      | sed s/XX_BALANCER_XX/$balancer/ \
      | sed s/XX_STIMBOX_XX/$stimBox/ \
      > $dirname/object.data
        
  cat tools/runMe.sh.proto \
      | sed s%XX_NNODES_XX%$nNodes% \
      | sed s%XX_NTASKS_XX%$nTasks% \
      | sed s%XX_MAXTIME_XX%$maxTime% \
      | sed s%XX_POOL_XX%$pool% \
      | sed s%XX_BANK_XX%$bank% \
      | sed s%XX_EXE_XX%$exe% \
      > $dirname/runMe.sh

done
done
done
