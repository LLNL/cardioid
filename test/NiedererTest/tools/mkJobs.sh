#!/bin/bash

# To Do

# 1. Might want to add a loop over different numbers of tasks to check
# that the results don't change with the decomposition.  

exe=../../../../src/cardioid-ibm_bgp_xlc 
pool=pdebug
maxTime=1:30
bank=dev
nNodes=64
nTasks=256


for dt in 0.05 0.01 0.005 0.001
do
for dx in 0.5 0.2 0.1 0.05
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
      nTasks=64
      nNodes=64
      xGrid=4; yGrid=4; zGrid=4
      tMax=150
      ;;
      0.2)
      stimBox=8
      nTasks=128
      nNodes=64
      xGrid=4; yGrid=4; zGrid=8
      tMax=60
      ;;
      0.1)
      stimBox=16
      nTasks=256
      nNodes=64
      xGrid=4; yGrid=8; zGrid=8
      tMax=50
      ;;
      0.05)
      stimBox=30
      nTasks=1024
      nNodes=256
      xGrid=8; yGrid=8; zGrid=16
      tMax=50
      ;;
      *)
      echo ERROR: undefined dx
      exit -1
  esac

  case $dt in
      0.05)
      snapshotRate=50
      ;;
      0.01)
      snapshotRate=250
      ;;
      0.005)
      snapshotRate=500
      ;;
      0.001)
      snapshotRate=2500
      ;;
      *)
      echo ERROR: undefined dt
      exit -1
  esac

  maxLoop=`echo $tMax/$dt | bc`

  cat tools/object.data.proto \
      | sed s/XX_DT_XX/$dt/ \
      | sed s/XX_DX_XX/$dx/ \
      | sed s/XX_MAXLOOP_XX/$maxLoop/ \
      | sed s/XX_SNAPSHOTRATE_XX/$snapshotRate/ \
      | sed s/XX_BALANCER_XX/$balancer/ \
      | sed s/XX_STIMBOX_XX/$stimBox/ \
      | sed s/XX_XGRID_XX/$xGrid/ \
      | sed s/XX_YGRID_XX/$yGrid/ \
      | sed s/XX_ZGRID_XX/$zGrid/ \
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
