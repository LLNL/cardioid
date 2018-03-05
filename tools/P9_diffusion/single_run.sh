export LD_LIBRARY_PATH=/home/ckim/project/CORAL/cardioid/serial_hpm/lib
#export OMP_DISPLAY_ENV=verbose
export OMP_NUM_THREADS=16
./diffusion -x 100 -y 62 -z 62 -K simd_cpu -N 200
