#!/bin/bash

# To Do


mpiExe=../../../..//bin/cardioid-bgq
spiExe=${mpiExe}-spi
spiTimeExe=${mpiExe}-spi-time
#pool=pdebug
#maxTime=1:30
#bank=dev



for halo in spitime
do
for jobSize in 1k 2k 4k 8k
do
for cellsPerNode in 150 300
do
for balancer in koradi grid
do

  dirname=run/${halo}_${cellsPerNode}_${balancer}_$jobSize
  if [ -d $dirname ] 
  then
      continue
  fi

  echo making $dirname
  mkdir -p $dirname

  case $halo in
      mpi)
      exe=$mpiExe
      ;;
      spi)
      exe=$spiExe
      ;;
      spitime)
      exe=$spiTimeExe
      ;;
  esac

  case $cellsPerNode in
      150)
      baseSize=10.7
      ;;
      300)
      baseSize=13.5
      ;;
      *)
      echo ERROR: undefined cellsPerNode
      exit -1
  esac



  case $jobSize in
      1k)
      nNodes=1024
      xGrid=8;   yGrid=8;   zGrid=16
      xSize=`echo $baseSize \* 1 | bc`
      ySize=`echo $baseSize \* 1 | bc`
      zSize=`echo $baseSize \* 2 | bc`
      ;;
      2k)
      nNodes=2048
      xGrid=8;   yGrid=16;  zGrid=16
      xSize=`echo $baseSize \* 1 | bc`
      ySize=`echo $baseSize \* 2 | bc`
      zSize=`echo $baseSize \* 2 | bc`
      ;;
      4k)
      nNodes=4096
      xGrid=16;  yGrid=16;  zGrid=16
      xSize=`echo $baseSize \* 2 | bc`
      ySize=`echo $baseSize \* 2 | bc`
      zSize=`echo $baseSize \* 2 | bc`
      ;;
      8k)
      nNodes=8192
      xGrid=16;  yGrid=16;  zGrid=32
      xSize=`echo $baseSize \* 2 | bc`
      ySize=`echo $baseSize \* 2 | bc`
      zSize=`echo $baseSize \* 4 | bc`
      ;;
      *)
      echo ERROR: undefined jobSize
      exit -1
  esac


#  maxLoop=`echo $tMax/$dt | bc`

  cat tools/object.data.proto \
      | sed s/XX_BALANCER_XX/$balancer/ \
      | sed s/XX_XGRID_XX/$xGrid/ \
      | sed s/XX_YGRID_XX/$yGrid/ \
      | sed s/XX_ZGRID_XX/$zGrid/ \
      | sed s/XX_XSIZE_XX/$xSize/ \
      | sed s/XX_YSIZE_XX/$ySize/ \
      | sed s/XX_ZSIZE_XX/$zSize/ \
      > $dirname/object.data
        
  nTasks=$nNodes
  cat tools/sbatchMe.sh.proto \
      | sed s%XX_NNODES_XX%$nNodes% \
      | sed s%XX_NTASKS_XX%$nTasks% \
      | sed s%XX_EXE_XX%$exe% \
      > $dirname/sbatchMe.sh

done
done
done
done
