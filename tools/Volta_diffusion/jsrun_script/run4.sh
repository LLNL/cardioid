#!/bin/bash 
nodes=1
ppn=4
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
#BSUB -R "affinity[core(11):distribute=pack]"
#BSUB -n ${nmpi}
#BSUB -x
#BSUB -q excl_ws_dd21
#BSUB -W 15
#---------------------------------------
ulimit -s 10240
ulimit -c 1000

/opt/ibm/spectrum_mpi/bin/mpirun -aff off -tag-output -np $nmpi ./snap 4rank.in mout.4rank.$tstamp

EOF
#---------------------------------------
bsub  <batch.job
