#!/bin/bash
#SBATCH --nodes=4096

export OMP_NUM_THREADS=64
export MUSPI_NUMBATIDS=203
export MUSPI_NUMINJFIFOS=3
export MUSPI_NUMRECFIFOS=3

srun --ntasks=4096 ./spi_test-bgq-spi
