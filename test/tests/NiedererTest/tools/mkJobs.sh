#!/bin/bash

# To Do

# 1. Might want to add a loop over different numbers of tasks to check
# that the results don't change with the decomposition.  

exe=../../../../bin/cardioid-bgq-spi
nNodes=64
nTasks=64
balancer=grid

for dt in 0.050 0.010 0.005 0.001
do
for dx in 0.50 0.20 0.10 0.05
do

  dirname=run/${balancer}_dt${dt}_dx${dx}
  if [ -d $dirname ] 
  then
      continue
  fi

  echo making $dirname
  mkdir -p $dirname

  case $dx in
      0.50)
      stimBox=4
      nTasks=64
      nNodes=64
      xGrid=4; yGrid=4; zGrid=4
      tMax=150
      ;;
      0.20)
      stimBox=8
      nTasks=64
      nNodes=64
      xGrid=4; yGrid=4; zGrid=4
      tMax=60
      ;;
      0.10)
      stimBox=16
      nTasks=64
      nNodes=64
      xGrid=4; yGrid=4; zGrid=4
      tMax=50
      ;;
      0.05)
      stimBox=30
      nTasks=64
      nNodes=64
      xGrid=4; yGrid=4; zGrid=4
      tMax=50
      ;;
      *)
      echo ERROR: undefined dx
      exit -1
  esac

  case $dt in
      0.050)
      snapshotRate=500000
      ;;
      0.010)
      snapshotRate=500000
      ;;
      0.005)
      snapshotRate=500000
      ;;
      0.001)
      snapshotRate=500000
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
        
  cat tools/sbatch.bgq.proto \
      | sed s%XX_NNODES_XX%$nNodes% \
      | sed s%XX_NTASKS_XX%$nTasks% \
      | sed s%XX_EXE_XX%$exe% \
      > $dirname/sbatch.bgq

done
done
