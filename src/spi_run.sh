#!/bin/bash -x
#export DCMF_EAGER=20000
#export DCMF_COLLECTIVES=0
#export DCMF_INTERRUPTS=1
#export DCMF_RECFIFO=bytes

#export BG_SHAREDMEMSIZE=64
export OMP_NUM_THREADS=64

# export PAMI_MAXCONTEXTS=1
# export PAMI_CONTEXT_POST=0
 export OMP_WAIT_POLICY=active
 export ATOMICS_OPT_LEVEL=5
# export LD_LIBRARY_PATH=/bgusr/hfwen/examples/NPB3.3-OMP/bin:/bgusr/hfwen/examples/hello:/gsa/rchgsa-p2/16/xlcmpbld/run/xlf/dev/bgq/daily/111124/opt/ibmcmp/lib64/bg:/bgsys/drivers/ppcfloor/comm/gcc/lib/:/bgsys/drivers/ppcfloor/comm/sys/lib
 export BLOCK=R01-M0-N02
# export TRACE_SEND_PATTERN=yes
# export PROFILE_BY_CALL_SITE=yes
# export SAVE_ALL_TASKS=yes
# export TRACE_ALL_EVENTS=yes
# export VPROF_PROFILE=yes
 export HPM_GROUP=2

export MUSPI_NUMBATIDS=203
export MUSPI_NUMINJFIFOS=3
export MUSPI_NUMRECFIFOS=3

#HPCT
#export TRACE_ALL_TASKS=yes
# export BLOCK=R00-M0-N00-64;R00-M0-N04-128;R00-M0-N10-64;R00-M0-N12-128;R00-M0-N08-64
#runjob --block $BLOCK --ranks-per-node 1 --np 1  --cwd $PWD --env_all : ./$1 >gtc.out 2>&1 &
runjob --label --block $BLOCK --ranks-per-node 1 --np 16  --cwd $PWD --env_all --exe ./$1 --args "object_FGR.data" > std.out_FGR
runjob --label --block $BLOCK --ranks-per-node 1 --np 16  --cwd $PWD --env_all --exe ./$1 --args "object_FLEX.data" > std.out_FLEX
grep '\[0\]' std.out_FGR > std.out_FGR_0
grep '\[0\]' std.out_FLEX > std.out_FLEX_0

