#!/bin/bash

#SBATCH --nodes=1024


#export BG_SHAREDMEMSIZE=64
#export PAMI_MAXCONTEXTS=1
#export PAMI_CONTEXT_POST=0
#export OMP_WAIT_POLICY=active
#export ATOMICS_OPT_LEVEL=5
#export HPM_GROUP=2

export OMP_NUM_THREADS=1
export MUSPI_NUMBATIDS=203
export MUSPI_NUMINJFIFOS=3
export MUSPI_NUMRECFIFOS=3

srun --ntasks=1024 /nfs/tmp2/emhm/src/Cardioid/EP/trunk/bin/cardioid-bgq-spi
