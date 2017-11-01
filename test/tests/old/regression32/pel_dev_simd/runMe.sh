#!/bin/sh

#PSUB -ln 6
#PSUB -tM 0:20
#PSUB -pool pbatch
#PSUB -b iscr

cd ${PSUB_SUBDIR:-.}

exe=../../../bin/cardioid-peloton

export OMP_NUM_THREADS=2

srun --ntasks=32 $exe

