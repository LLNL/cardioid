#!/bin/bash 
nodes=1
ppn=6
let nmpi=$nodes*$ppn
tstamp=`date +%m_%d_%H_%M_%S`
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=${ppn}]"
#BSUB -R "select[ngpus=6] rusage[ngpus_shared=20]"
#BSUB -R "select[type=RHEL7_4]" 
#BSUB -env "LSB_START_JOB_MPS=N"
#BSUB -R "affinity[core(2,same=socket):cpubind=socket:distribute=balance]"
#BSUB -n ${nmpi}
#BSUB -x
#BSUB -q excl_ws_dd10
##BSUB -q excl_ws_dd10_int
##BSUB -q shared_ws_dd10
##BSUB -q shared_ws_dd10_int
#BSUB -W 15
#---------------------------------------
ulimit -s 10240
ulimit -c 1000


#more /proc/cpuinfo
/home/walkup/openmpi-2.1.1/gnu-4.8.5/bin/mpirun --mca orte_base_help_aggregate 0 -np $nmpi -tag-output -aff off ./snap 6rank.in mout.6rank.$tstamp


#export LD_PRELOAD=/home/walkup/mpitrace/spectrum_mpi/libmpitrace.so
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi /usr/local/cuda/bin/nvprof --log-file nvN1.out -m inst_per_warp ./snap 1rank.in mout.1rank.$tstamp
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi /usr/local/cuda/bin/nvprof --print-gpu-trace --track-memory-allocations on --log-file nvN1.out ./snap 1rank.in mout.1rank.$tstamp
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi /usr/local/cuda/bin/nvprof --log-file nvN1.out ./snap 1rank.in mout.1rank.$tstamp
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi valgrind --tool=memcheck ./snap 1rank.in mout.1rank.$tstamp
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi /usr/local/cuda/bin/cuda-memcheck ./snap 1rank.in mout.1rank.$tstamp
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi gdb ./snap
#unset LD_PRELOAD
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi -tag-output ./snap 1rank.in mout.1rank.$tstamp
#/home/walkup/openmpi-2.1.1/gnu-4.8.5/bin/mpirun -np $nmpi -tag-output ./snap 1rank.in mout.1rank.$tstamp
#/home/walkup/openmpi-2.1.1/gnu-4.8.5/bin/mpirun -np $nmpi ./snap 1rank.in mout.1rank.$tstamp
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi /usr/local/cuda/bin/cuda-memcheck --tool racecheck ./snap 1rank.in mout.1rank.$tstamp
EOF
#---------------------------------------
bsub  <batch.job
