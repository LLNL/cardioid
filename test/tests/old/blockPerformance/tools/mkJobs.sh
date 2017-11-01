#!/bin/bash

mpiExe=../../../..//bin/cardioid-bgq
spiExe=${mpiExe}-spi
ossExe=${mpiExe}-spi-oss
hpmExe=${mpiExe}-spi-hpm


for jobType in spi hpm
do
for jobSize in 1  32 1k 2k 4k 8k 16k
do
for cellsPerNode in 224
do
for balancer in grid
do

  dirname=run/${jobType}_${cellsPerNode}_${balancer}_$jobSize
  if [ -d $dirname ] 
  then
      continue
  fi
  ossdir=`pwd`/$dirname/rawoss_data

  echo making $dirname
  mkdir -p $dirname
  mkdir -p $ossdir

  case $jobType in
      mpi)
      exe=$mpiExe
      ;;
      spi)
      exe=$spiExe
      ;;
      hpm)
      exe=$hpmExe
      ;;
      oss)
      exe=$ossExe
      ;;
  esac

  case $cellsPerNode in
      150)
      xbaseSize=15
      ybaseSize=16
      zbaseSize=12.5
      ;;
      224)
      xbaseSize=16
      ybaseSize=16
      zbaseSize=14
      ;;
      300)
      xbaseSize=20
      ybaseSize=16
      zbaseSize=15
      ;;
      *)
      echo ERROR: undefined cellsPerNode
      exit -1
  esac



  case $jobSize in
      1)
      nNodes=
      xGrid=1;   yGrid=1;   zGrid=1
      nNodeName=$jobSize
      ;;
      32)
      nNodes=
      xGrid=2;   yGrid=4;   zGrid=4
      nNodeName=$jobSize
      ;;
      1k)
      nNodes=1024
      xGrid=8;   yGrid=8;   zGrid=16
      nNodeName=$jobSize
      ;;
      2k)
      nNodes=2048
      xGrid=16;   yGrid=8;  zGrid=16
      nNodeName=$jobSize
      ;;
      4k)
      nNodes=4096
      xGrid=16;  yGrid=16;  zGrid=16
      nNodeName=$jobSize
      ;;
      8k)
      nNodes=8192
      xGrid=16;  yGrid=16;  zGrid=32
      nNodeName=$jobSize
      ;;
      16k)
      nNodes=16384
      xGrid=32;  yGrid=16;  zGrid=32
      nNodeName=$jobSize
      ;;
      24k)
      nNodes=24576
      xGrid=32;  yGrid=24;  zGrid=32
      nNodeName=$jobSize
      ;;
      48k)
      nNodes=49152
      xGrid=32;  yGrid=48;  zGrid=32
      nNodeName=$jobSize
      ;;
      *)
      echo ERROR: undefined jobSize
      exit -1
  esac
  xSize=`echo $xbaseSize \* $xGrid | bc`
  ySize=`echo $ybaseSize \* $yGrid | bc`
  zSize=`echo $zbaseSize \* $zGrid | bc`


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
      | sed s%XX_OSSDIR_XX%$ossdir% \
      | sed s%XX_EXE_XX%$exe% \
      > $dirname/sbatchMe.sh

done
done
done
done
