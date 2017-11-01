#!/bin/sh

#PSUB -ln 3
#PSUB -tM 0:20
#PSUB -pool pbatch
#PSUB -b iscr

cd ${PSUB_SUBDIR:-.}

exe=../../../bin/cardioid-peloton

export OMP_NUM_THREADS=1

srun --ntasks=32 $exe

