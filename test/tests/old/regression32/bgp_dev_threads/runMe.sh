#!/bin/sh

#PSUB -ln 64
#PSUB -tM 1:30
#PSUB -pool pdebug
#PSUB -b dev

cd ${PSUB_SUBDIR:-.}

exe=../../../bin/cardioid-bgp

mpirun -mode SMP -np 32 -exe $exe -env OMP_NUM_THREADS=4
