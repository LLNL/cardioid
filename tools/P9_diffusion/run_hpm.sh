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
#BSUB -n ${nmpi}
#BSUB -x
#BSUB -q excl
#BSUB -W 15
#---------------------------------------
ulimit -s 10240
export LD_LIBRARY_PATH=/home/ckim/project/CORAL/cardioid/serial_hpm/lib
export OMP_PLACES={0:16:4}
export OMP_WAIT_POLICY=active
./diffusion -x 100 -y 62 -z 62 -K simd_cpu_thr2 -N 100
./diffusion -x 100 -y 62 -z 62 -K simd_cpu_thr2 -N 100

./diffusion -x 200 -y 62 -z 62 -K simd_cpu_thr2 -N 100
./diffusion -x 200 -y 62 -z 62 -K simd_cpu_thr2 -N 100

export OMP_PLACES={0:8:8}
export OMP_WAIT_POLICY=active
./diffusion -x 200 -y 32 -z 62 -K simd_cpu_thr2 -N 100
./diffusion -x 200 -y 32 -z 62 -K simd_cpu_thr2 -N 100

export OMP_PLACES={0:9:8}
export OMP_WAIT_POLICY=active
./diffusion -x 200 -y 47 -z 47 -K simd_cpu_thr2 -N 100
./diffusion -x 200 -y 47 -z 47 -K simd_cpu_thr2 -N 100
EOF
#---------------------------------------
bsub  <batch.job
