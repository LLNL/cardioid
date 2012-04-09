#!/bin/sh

#PSUB -ln 64
#PSUB -tM 0:30
#PSUB -pool pdebug
#PSUB -b science
#SBATCH -N 32

cd ${PSUB_SUBDIR:-.}

exe=../../../bin/cardioid-bgq-spi

export OMP_NUM_THREADS=64
export MUSPI_NUMBATIDS=203
export MUSPI_NUMINJFIFOS=3
export MUSPI_NUMRECFIFOS=3

srun -N 32 -n 32 $exe

