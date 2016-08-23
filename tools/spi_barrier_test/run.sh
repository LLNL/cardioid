#!/bin/bash -x
#export BG_SHAREDMEMSIZE=48
#export PAMI_MAXCONTEXTS=2
#export PAMI_CONTEXT_POST=1
#export PAMI_COMM_THREADS=0
export MUSPI_NUMBATIDS=203
export MUSPI_NUMINJFIFOS=3
export MUSPI_NUMRECFIFOS=3
export BLOCK=R03-M0
runjob --label --block ${BLOCK} --cwd $PWD --ranks-per-node 1  --np 512 --env-all : spi_test-spi  > stdout
