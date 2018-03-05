#nvcc -arch=sm_60 -Dflux_dump -c --maxrregcount=63 CUDADiffusion_ker.cu
nvcc -Xptxas -v -arch=sm_70 -c --maxrregcount=127 CUDADiffusion_ker.cu
nvcc -arch=sm_70 main.cu CUDADiffusion_ker.o -o diff.x
#nvprof ./diff.x
#nvprof -m flop_count_dp ./diff.x
