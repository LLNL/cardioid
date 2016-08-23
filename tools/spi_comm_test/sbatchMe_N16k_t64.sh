#!/bin/bash
#SBATCH --nodes=16384

export OMP_NUM_THREADS=64
export MUSPI_NUMBATIDS=203
export MUSPI_NUMINJFIFOS=3
export MUSPI_NUMRECFIFOS=3

srun --ntasks=16384 ./spi_test-bgq-spi
