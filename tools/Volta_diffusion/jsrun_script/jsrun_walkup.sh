#!/bin/bash 
nodes=10
ppn=6
let nmpi=$nodes*$ppn
#--------------------------------------
cat >batch.job <<EOF
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -nnodes ${nodes}
#BSUB -q excl 
#BSUB -W 15
#---------------------------------
ulimit -s 10240
export OMP_NUM_THREADS=1
export BIND_SLOTS=24
export RANKS_PER_NODE=${ppn}
export TRACE_ALL_EVENTS=yes
export TRACE_ALL_TASKS=yes
export LD_PRELOAD=/shared/ext_sw/tools/mpitrace/spectrum_mpi/libmpitrace.so
 /opt/ibm/spectrum_mpi/jsm_pmix/bin/jsrun --smpiargs -gpu --rs_per_host 1 --tasks_per_rs ${ppn} \
 --cpu_per_rs 44 --gpu_per_rs 6 --nrs ${nodes} -d plane:${ppn} ./helper_lammps_walkup.sh  \
 ../src/lmp_kokkos_cuda_mpi -k on g 6 -sf kk -pk kokkos \
 neigh half neigh/qeq full newton on -v x 60 -v y 48 -v z 36 -in in.reaxc.hns_543 -nocite

EOF
#---------------------------------------
 bsub <batch.job
