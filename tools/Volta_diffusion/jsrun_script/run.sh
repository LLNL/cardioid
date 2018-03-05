#!/bin/bash 
nodes=1
ppn=1
let nmpi=$nodes*$ppn
tstamp=`date +%m_%d_%H_%M_%S`
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=${ppn}]"
#BSUB -R "select[ngpus=4] rusage[ngpus_shared=20]"
#BSUB -env "LSB_START_JOB_MPS=N"
#BSUB -n ${nmpi}
#BSUB -x
#BSUB -q excl
#BSUB -W 15
#---------------------------------------
ulimit -s 10240
ulimit -c 1000
#export LD_PRELOAD=/home/walkup/mpitrace/spectrum_mpi/libmpitrace.so
#/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi /usr/local/cuda/bin/nvprof --log-file nv1.out -m inst_per_warp ./snap 1rank.in mout.1rank.$tstamp
/opt/ibm/spectrum_mpi/bin/mpirun -aff off -np $nmpi /usr/local/cuda/bin/nvprof --log-file nv1.out ./snap 1rank.in mout.1rank.$tstamp
#unset LD_PRELOAD
EOF
#---------------------------------------
bsub  <batch.job
