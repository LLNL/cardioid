#!/bin/bash 
nodes=1
ppn=2
let nmpi=$nodes*$ppn
tstamp=`date +%m_%d_%H_%M_%S`
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "span[ptile=${ppn}]"
#BSUB -R "select[ngpus=4] rusage[ngpus_shared=20]"
#BSUB -env "LSB_START_JOB_MPS=N"
#BSUB -R "affinity[core(5):distribute=pack]"
#BSUB -n ${nmpi}
#BSUB -x
#BSUB -q excl
#BSUB -W 15
#---------------------------------------
ulimit -s 10240
ulimit -c 1000
export LD_PRELOAD=/home/walkup/mpitrace/spectrum_mpi/libmpitrace.so
/opt/ibm/spectrum_mpi/bin/mpirun -tag-output -aff off -np $nmpi /usr/local/cuda/bin/nvprof --log-file nv%q{OMPI_COMM_WORLD_RANK}.out ./snap 2rank.in mout.2rank.$tstamp
unset LD_PRELOAD
EOF
#---------------------------------------
bsub  <batch.job
