#!/bin/bash

# To Do

# 1. Might want to add a loop over different numbers of tasks to check
# that the results don't change with the decomposition.  

exe=../../../../bin/cardioid-bgq-spi
nNodes=64
nTasks=64
balancer=grid

for dt in 0.05 0.01 0.005
do
for dx in 0.25 0.10 0.05
do

  dirname=run/${balancer}_dt${dt}_dx${dx}
  if [ -d $dirname ] 
  then
      continue
  fi

  echo making $dirname
  mkdir -p $dirname

  case $dx in
      0.25)
      s1Box=6
      s2BoxX=5
      s2BoxY=21
      s2BoxZ=73 
      nTasks=64
      nNodes=64
      xGrid=4; yGrid=4; zGrid=4
      tMax=1000
      ;;
      0.10)
      s1Box=15
      s2BoxX=14
      s2BoxY=54
      s2BoxZ=184 
      nTasks=64
      nNodes=64
      xGrid=4; yGrid=4; zGrid=4
      tMax=1000
      ;;
      0.05)
      s1Box=30
      s2BoxX=29
      s2BoxY=109
      s2BoxZ=369 
      nTasks=64
      nNodes=64
      xGrid=4; yGrid=4; zGrid=4
      tMax=1000
      ;;
      *)
      echo ERROR: undefined dx
      exit -1
  esac

  maxLoop=`echo $tMax/$dt | bc`

  cat tools/object.data.proto \
      | sed s/XX_DT_XX/$dt/ \
      | sed s/XX_DX_XX/$dx/ \
      | sed s/XX_MAXLOOP_XX/$maxLoop/ \
      | sed s/XX_BALANCER_XX/$balancer/ \
      | sed s/XX_S1BOX_XX/$s1Box/ \
      | sed s/XX_S2BOX_X_XX/$s2BoxX/ \
      | sed s/XX_S2BOX_Y_XX/$s2BoxY/ \
      | sed s/XX_S2BOX_Z_XX/$s2BoxZ/ \
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
