nvcc -arch=sm_60 -Dflux_dump -c --maxrregcount=63 CUDADiffusion_ker.cu
nvcc -arch=sm_60 main.cu CUDADiffusion_ker.o -o diff.x
#nvprof ./diff.x
#nvprof -m flop_count_dp ./diff.x
