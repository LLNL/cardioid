nvcc -arch=sm_60 CUDADiffusion_ker.cu -o diff.x
nvprof ./diff.x
nvprof -m flop_count_dp ./diff.x
